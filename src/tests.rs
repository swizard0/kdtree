use std::cmp::{min, max, Ordering};
use rand::{self, Rng};
use super::{KdvTree, Intersection, NearestShape};

#[derive(Clone, Copy, PartialEq, Debug)]
struct Point2d { x: i32, y: i32, }

#[derive(Clone, Debug)]
enum Axis { X, Y, }

fn cmp_points(axis: &Axis, a: &Point2d, b: &Point2d) -> Ordering {
    match axis {
        &Axis::X => a.x.cmp(&b.x),
        &Axis::Y => a.y.cmp(&b.y),
    }
}

#[derive(PartialEq, Clone, Debug)]
struct Rect2d { lt: Point2d, rb: Point2d, }

impl super::BoundingVolume<Point2d> for Rect2d {
    fn min_corner(&self) -> Point2d { self.lt }
    fn max_corner(&self) -> Point2d { self.rb }
}

#[derive(PartialEq, Debug)]
struct Line2d { src: Point2d, dst: Point2d, }

fn get_bounding_volume(shape: &Line2d) -> Rect2d {
    Rect2d {
        lt: Point2d { x: min(shape.src.x, shape.dst.x), y: min(shape.src.y, shape.dst.y), },
        rb: Point2d { x: max(shape.src.x, shape.dst.x), y: max(shape.src.y, shape.dst.y), },
    }
}

fn get_cut_point(cut_axis: &Axis, points: &mut Iterator<Item = Point2d>) -> Option<Point2d> {
    if let Some(first_point) = points.next() {
        let mut point_min = Point2d { x: first_point.x, y: first_point.y, };
        let mut point_max = Point2d { x: first_point.x, y: first_point.y, };
        let mut total_x = first_point.x as i64;
        let mut total_y = first_point.y as i64;
        let mut total = 1;
        for p in points {
            total_x += p.x as i64;
            total_y += p.y as i64;
            total += 1;
            if p.x < point_min.x { point_min.x = p.x; }
            if p.y < point_min.y { point_min.y = p.y; }
            if p.x > point_max.x { point_max.x = p.x; }
            if p.y > point_max.y { point_max.y = p.y; }
        }
        let point_mid = Point2d {
            x: (total_x as f64 / total as f64) as i32,
            y: (total_y as f64 / total as f64) as i32,
        };
        match cut_axis {
            &Axis::X if point_mid.x == point_min.x || point_mid.x == point_max.x =>
                None,
            &Axis::Y if point_mid.y == point_min.y || point_mid.y == point_max.y =>
                None,
            _ =>
                Some(point_mid)
        }
    } else {
        None
    }
}

struct Cutter {
    cut_limit: i32,
}

impl super::BoundingVolumesCutter<Axis, Point2d, Rect2d, Line2d> for Cutter {
    type Error = ();

    fn cut(&mut self, shape: &Line2d, fragment: &Rect2d, cut_axis: &Axis, cut_point: &Point2d) ->
        Result<Option<(Rect2d, Rect2d)>, ()>
    {
        match cut_axis {
            &Axis::X => if cut_point.x >= fragment.lt.x && cut_point.x <= fragment.rb.x {
                if fragment.rb.x - fragment.lt.x < self.cut_limit {
                    Ok(None)
                } else {
                    let factor = (cut_point.x - shape.src.x) as f64 / (shape.dst.x - shape.src.x) as f64;
                    let y = shape.src.y as f64 + (factor * (shape.dst.y - shape.src.y) as f64);
                    let left_point = if shape.src.x < shape.dst.x { shape.src } else { shape.dst };
                    let left_bound = Rect2d {
                        lt: Point2d {
                            x: fragment.lt.x,
                            y: if left_point.y < y as i32 { fragment.lt.y } else { y as i32 },
                        },
                        rb: Point2d {
                            x: cut_point.x,
                            y: if left_point.y < y as i32 { y as i32 } else { fragment.rb.y },
                        }
                    };
                    let right_point = if shape.src.x < shape.dst.x { shape.dst } else { shape.src };
                    let right_bound = Rect2d {
                        lt: Point2d {
                            x: cut_point.x,
                            y: if right_point.y < y as i32 { fragment.lt.y } else { y as i32 },
                        },
                        rb: Point2d {
                            x: fragment.rb.x,
                            y: if right_point.y < y as i32 { y as i32 } else { fragment.rb.y },
                        },
                    };
                    Ok(Some((left_bound, right_bound)))
                }
            } else {
                return Ok(None);
            },
            &Axis::Y => if cut_point.y >= fragment.lt.y && cut_point.y <= fragment.rb.y {
                if fragment.rb.y - fragment.lt.y < self.cut_limit {
                    Ok(None)
                } else {
                    let factor = (cut_point.y - shape.src.y) as f64 / (shape.dst.y - shape.src.y) as f64;
                    let x = shape.src.x as f64 + (factor * (shape.dst.x - shape.src.x) as f64);
                    let upper_point = if shape.src.y < shape.dst.y { shape.src } else { shape.dst };
                    let upper_bound = Rect2d {
                        lt: Point2d {
                            x: if upper_point.x < x as i32 { fragment.lt.x } else { x as i32 },
                            y: fragment.lt.y,
                        },
                        rb: Point2d {
                            x: if upper_point.x < x as i32 { x as i32 } else { fragment.rb.x },
                            y: cut_point.y,
                        }
                    };
                    let lower_point = if shape.src.y < shape.dst.y { shape.dst } else { shape.src };
                    let lower_bound = Rect2d {
                        lt: Point2d {
                            x: if lower_point.x < x as i32 { fragment.lt.x } else { x as i32 },
                            y: cut_point.y,
                        },
                        rb: Point2d {
                            x: if lower_point.x < x as i32 { x as i32 } else { fragment.rb.x },
                            y: fragment.rb.y,
                        },
                    };
                    Ok(Some((upper_bound, lower_bound)))
                }
            } else {
                return Ok(None);
            },
        }
    }
}

fn bv_to_cut_point_sq_dist(axis: &Axis, bounding_volume: &Rect2d, cut_point: &Point2d) -> i64 {
    let dist = match axis {
        &Axis::X =>
            min((bounding_volume.lt.x - cut_point.x).abs() as i64, (bounding_volume.rb.x - cut_point.x).abs() as i64),
        &Axis::Y =>
            min((bounding_volume.lt.y - cut_point.y).abs() as i64, (bounding_volume.rb.y - cut_point.y).abs() as i64),
    };
    dist * dist
}

fn bv_to_bv_sq_dist(bounding_volume_a: &Rect2d, bounding_volume_b: &Rect2d) -> i64 {
    let center_a = Point2d {
        x: (bounding_volume_a.lt.x + bounding_volume_a.rb.x) / 2,
        y: (bounding_volume_a.lt.y + bounding_volume_a.rb.y) / 2,
    };
    let center_b = Point2d {
        x: (bounding_volume_b.lt.x + bounding_volume_b.rb.x) / 2,
        y: (bounding_volume_b.lt.y + bounding_volume_b.rb.y) / 2,
    };
    ((center_b.x - center_a.x) as i64 * (center_b.x - center_a.x) as i64) + ((center_b.y - center_a.y) as i64 * (center_b.y - center_a.y) as i64)
}

fn random_line<R>(rng: &mut R) -> Line2d where R: Rng {
    let src = Point2d { x: rng.gen_range(-1000000, 1000000), y: rng.gen_range(-1000000, 1000000), };
    let dst = Point2d { x: src.x + rng.gen_range(-100000, 100000), y: src.y + rng.gen_range(-100000, 100000), };
    Line2d { src, dst }
}

#[test]
fn kdv_tree_basic() {
    let shapes = vec![Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, }];
    let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, cmp_points, get_bounding_volume, get_cut_point, Cutter { cut_limit: 10, }).unwrap();

    assert_eq!(
        tree.intersects(
            &Line2d { src: Point2d { x: 116, y: 116, }, dst: Point2d { x: 180, y: 180, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
            .collect::<Result<Vec<_>, _>>(),
        Ok(vec![])
    );
    assert_eq!(
        tree.intersects(
            &Line2d { src: Point2d { x: 32, y: 48, }, dst: Point2d { x: 48, y: 64, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
            .collect::<Result<Vec<_>, _>>(),
        Ok(vec![])
    );
    assert_eq!(
        tree.intersects(
            &Line2d { src: Point2d { x: 48, y: 32, }, dst: Point2d { x: 64, y: 48, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
            .collect::<Result<Vec<_>, _>>(),
        Ok(vec![])
    );

    let intersects: Result<Vec<_>, _> = tree
        .intersects(
            &Line2d { src: Point2d { x: 16, y: 64, }, dst: Point2d { x: 80, y: 64, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
        .collect();
    assert_eq!(intersects, Ok(vec![
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 64, y: 64 }, rb: Point2d { x: 72, y: 72 } },
            needle_fragment: Rect2d { lt: Point2d { x: 66, y: 64 }, rb: Point2d { x: 72, y: 64 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 64, y: 64 }, rb: Point2d { x: 72, y: 72 } },
            needle_fragment: Rect2d { lt: Point2d { x: 60, y: 64 }, rb: Point2d { x: 66, y: 64 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 56, y: 56 }, rb: Point2d { x: 64, y: 64 } },
            needle_fragment: Rect2d { lt: Point2d { x: 62, y: 64 }, rb: Point2d { x: 68, y: 64 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 56, y: 56 }, rb: Point2d { x: 64, y: 64 } },
            needle_fragment: Rect2d { lt: Point2d { x: 56, y: 64 }, rb: Point2d { x: 62, y: 64 } },
        },
    ]));
    let intersects: Result<Vec<_>, _> = tree
        .intersects(
            &Line2d { src: Point2d { x: 64, y: 16, }, dst: Point2d { x: 64, y: 80, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
        .collect();
    assert_eq!(intersects, Ok(vec![
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 64, y: 64 }, rb: Point2d { x: 72, y: 72 } },
            needle_fragment: Rect2d { lt: Point2d { x: 64, y: 72 }, rb: Point2d { x: 64, y: 80 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 64, y: 64 }, rb: Point2d { x: 72, y: 72 } },
            needle_fragment: Rect2d { lt: Point2d { x: 64, y: 64 }, rb: Point2d { x: 64, y: 72 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 56, y: 56 }, rb: Point2d { x: 64, y: 64 } },
            needle_fragment: Rect2d { lt: Point2d { x: 64, y: 58 }, rb: Point2d { x: 64, y: 64 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 56, y: 56 }, rb: Point2d { x: 64, y: 64 } },
            needle_fragment: Rect2d { lt: Point2d { x: 64, y: 52 }, rb: Point2d { x: 64, y: 58 } },
        },
    ]));
}

#[test]
fn kdv_tree_triangle() {
    let shapes = vec![
        Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 16, }, },
        Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, },
        Line2d { src: Point2d { x: 80, y: 16, }, dst: Point2d { x: 80, y: 80, }, },
    ];
    let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, cmp_points, get_bounding_volume, get_cut_point, Cutter { cut_limit: 10, }).unwrap();

    assert_eq!(
        tree.intersects(
            &Line2d { src: Point2d { x: 70, y: 45, }, dst: Point2d { x: 75, y: 50, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
            .collect::<Result<Vec<_>, _>>(),
        Ok(vec![])
    );

    let intersects: Result<Vec<_>, _> = tree
        .intersects(
            &Line2d { src: Point2d { x: 8, y: 48, }, dst: Point2d { x: 88, y: 48, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
        .collect();
    assert_eq!(intersects, Ok(vec![
        Intersection {
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 44 }, rb: Point2d { x: 80, y: 69 } },
            needle_fragment: Rect2d { lt: Point2d { x: 74, y: 48 }, rb: Point2d { x: 81, y: 48 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 42, y: 42 }, rb: Point2d { x: 50, y: 50 } },
            needle_fragment: Rect2d { lt: Point2d { x: 50, y: 48 }, rb: Point2d { x: 58, y: 48 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 42, y: 42 }, rb: Point2d { x: 50, y: 50 } },
            needle_fragment: Rect2d { lt: Point2d { x: 42, y: 48 }, rb: Point2d { x: 50, y: 48 } },
        },
    ]));
    let intersects: Result<Vec<_>, _> = tree
        .intersects(
            &Line2d { src: Point2d { x: 40, y: 10, }, dst: Point2d { x: 90, y: 60, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 10, },
        )
        .collect();
    assert_eq!(intersects, Ok(vec![
        Intersection {
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 44 }, rb: Point2d { x: 80, y: 69 } },
            needle_fragment: Rect2d { lt: Point2d { x: 74, y: 44 }, rb: Point2d { x: 82, y: 52 } },
        },
        Intersection {
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 29, y: 16 }, rb: Point2d { x: 58, y: 16 } },
            needle_fragment: Rect2d { lt: Point2d { x: 40, y: 10 }, rb: Point2d { x: 48, y: 18 } },
        },
    ]));
}

#[test]
fn kdv_tree_nearest() {
    let shapes = vec![
        Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 16, }, },
        Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, },
        Line2d { src: Point2d { x: 80, y: 16, }, dst: Point2d { x: 80, y: 80, }, },
    ];
    let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, cmp_points, get_bounding_volume, get_cut_point, Cutter { cut_limit: 10, }).unwrap();

    let nearest: Result<Vec<_>, _> = tree
        .nearest(
            &Line2d { src: Point2d { x: 44, y: 8, }, dst: Point2d { x: 52, y: 12, }, },
            cmp_points,
            get_bounding_volume,
            Cutter { cut_limit: 10, },
            bv_to_cut_point_sq_dist,
            bv_to_bv_sq_dist,
        )
        .collect();
    assert_eq!(nearest, Ok(vec![
        NearestShape {
            dist: 61,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 29, y: 16 }, rb: Point2d { x: 58, y: 16 } },
        },
        NearestShape {
            dist: 360,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 58, y: 16 }, rb: Point2d { x: 74, y: 16 } },
        },
        NearestShape {
            dist: 612,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 19, y: 16 }, rb: Point2d { x: 29, y: 16 } },
        },
        NearestShape {
            dist: 724,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 26, y: 26 }, rb: Point2d { x: 34, y: 34 } },
        },
        NearestShape {
            dist: 820,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 18, y: 18 }, rb: Point2d { x: 26, y: 26 } },
        },
        NearestShape {
            dist: 877,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 74, y: 16 }, rb: Point2d { x: 80, y: 16 } },
        },
        NearestShape {
            dist: 884,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 34, y: 34 }, rb: Point2d { x: 42, y: 42 } },
        },
        NearestShape {
            dist: 997,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 16, y: 16 }, rb: Point2d { x: 19, y: 16 } },
        },
        NearestShape {
            dist: 1010,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 16, y: 16 }, rb: Point2d { x: 18, y: 18 } },
        },
        NearestShape {
            dist: 1300,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 42, y: 42 }, rb: Point2d { x: 50, y: 50 } },
        },
        NearestShape {
            dist: 1553,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 23 }, rb: Point2d { x: 80, y: 44 } },
        },
        NearestShape {
            dist: 1972,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 50, y: 50 }, rb: Point2d { x: 58, y: 58 } },
        },
        NearestShape {
            dist: 2900,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 58, y: 58 }, rb: Point2d { x: 66, y: 66 } },
        },
        NearestShape {
            dist: 3140,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 44 }, rb: Point2d { x: 80, y: 69 } },
        },
        NearestShape {
            dist: 4084,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 66, y: 66 }, rb: Point2d { x: 74, y: 74 } },
        },
        NearestShape {
            dist: 5330,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 74, y: 74 }, rb: Point2d { x: 80, y: 80 } },
        },
    ]));
    let nearest: Result<Vec<_>, _> = tree
        .nearest(
            &Line2d { src: Point2d { x: 44, y: 88, }, dst: Point2d { x: 52, y: 92, }, },
            cmp_points,
            get_bounding_volume,
            Cutter { cut_limit: 10, },
            bv_to_cut_point_sq_dist,
            bv_to_bv_sq_dist,
        )
        .collect();
    assert_eq!(nearest, Ok(vec![
        NearestShape {
            dist: 884,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 66, y: 66 }, rb: Point2d { x: 74, y: 74 } },
        },
        NearestShape {
            dist: 980,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 58, y: 58 }, rb: Point2d { x: 66, y: 66 } },
        },
        NearestShape {
            dist: 1010,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 74, y: 74 }, rb: Point2d { x: 80, y: 80 } },
        },
        NearestShape {
            dist: 1332,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 50, y: 50 }, rb: Point2d { x: 58, y: 58 } },
        },
        NearestShape {
            dist: 1940,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 42, y: 42 }, rb: Point2d { x: 50, y: 50 } },
        },
        NearestShape {
            dist: 2180,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 44 }, rb: Point2d { x: 80, y: 69 } },
        },
        NearestShape {
            dist: 2804,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 34, y: 34 }, rb: Point2d { x: 42, y: 42 } },
        },
        NearestShape {
            dist: 3924,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 26, y: 26 }, rb: Point2d { x: 34, y: 34 } },
        },
        NearestShape {
            dist: 4273,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 23 }, rb: Point2d { x: 80, y: 44 } },
        },
        NearestShape {
            dist: 5300,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 18, y: 18 }, rb: Point2d { x: 26, y: 26 } },
        },
        NearestShape {
            dist: 5501,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 29, y: 16 }, rb: Point2d { x: 58, y: 16 } },
        },
        NearestShape {
            dist: 5800,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 58, y: 16 }, rb: Point2d { x: 74, y: 16 } },
        },
        NearestShape {
            dist: 6052,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 19, y: 16 }, rb: Point2d { x: 29, y: 16 } },
        },
        NearestShape {
            dist: 6290,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 16, y: 16 }, rb: Point2d { x: 18, y: 18 } },
        },
        NearestShape {
            dist: 6317,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 74, y: 16 }, rb: Point2d { x: 80, y: 16 } },
        },
        NearestShape {
            dist: 6437,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 16, y: 16 }, rb: Point2d { x: 19, y: 16 } },
        },
    ]));
    let nearest: Result<Vec<_>, _> = tree
        .nearest(
            &Line2d { src: Point2d { x: 144, y: 48, }, dst: Point2d { x: 152, y: 52, }, },
            cmp_points,
            get_bounding_volume,
            Cutter { cut_limit: 10, },
            bv_to_cut_point_sq_dist,
            bv_to_bv_sq_dist,
        )
        .collect();
    assert_eq!(nearest, Ok(vec![
        NearestShape {
            dist: 4660,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 44 }, rb: Point2d { x: 80, y: 69 } },
        },
        NearestShape {
            dist: 4913,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 23 }, rb: Point2d { x: 80, y: 44 } },
        },
        NearestShape {
            dist: 5065,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 69 }, rb: Point2d { x: 80, y: 74 } },
        },
        NearestShape {
            dist: 5353,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 74 }, rb: Point2d { x: 80, y: 80 } },
        },
        NearestShape {
            dist: 5585,
            shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 80, y: 16 }, rb: Point2d { x: 80, y: 23 } },
        },
        NearestShape {
            dist: 5770,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 74, y: 74 }, rb: Point2d { x: 80, y: 80 } },
        },
        NearestShape {
            dist: 6197,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 74, y: 16 }, rb: Point2d { x: 80, y: 16 } },
        },
        NearestShape {
            dist: 6484,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 66, y: 66 }, rb: Point2d { x: 74, y: 74 } },
        },
        NearestShape {
            dist: 7540,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 58, y: 58 }, rb: Point2d { x: 66, y: 66 } },
        },
        NearestShape {
            dist: 7880,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 58, y: 16 }, rb: Point2d { x: 74, y: 16 } },
        },
        NearestShape {
            dist: 8852,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 50, y: 50 }, rb: Point2d { x: 58, y: 58 } },
        },
        NearestShape {
            dist: 10420,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 42, y: 42 }, rb: Point2d { x: 50, y: 50 } },
        },
        NearestShape {
            dist: 12181,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 29, y: 16 }, rb: Point2d { x: 58, y: 16 } },
        },
        NearestShape {
            dist: 12244,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 34, y: 34 }, rb: Point2d { x: 42, y: 42 } },
        },
        NearestShape {
            dist: 14324,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 26, y: 26 }, rb: Point2d { x: 34, y: 34 } },
        },
        NearestShape {
            dist: 16532,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 19, y: 16 }, rb: Point2d { x: 29, y: 16 } },
        },
        NearestShape {
            dist: 16660,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 18, y: 18 }, rb: Point2d { x: 26, y: 26 } },
        },
        NearestShape {
            dist: 18250,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 16, y: 16 }, rb: Point2d { x: 18, y: 18 } },
        },
        NearestShape {
            dist: 18317,
            shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } },
            shape_fragment: &Rect2d { lt: Point2d { x: 16, y: 16 }, rb: Point2d { x: 19, y: 16 } },
        },
    ]));
}

#[test]
fn kdv_tree_nearest_check() {
    let tree = KdvTree::build(
        vec![Axis::X, Axis::Y],
        (0 .. 1000).map(|_| random_line(&mut rand::thread_rng())),
        cmp_points,
        get_bounding_volume,
        get_cut_point,
        Cutter { cut_limit: 1000, }
    ).unwrap();

    let mut sample = random_line(&mut rand::thread_rng());
    sample.dst.x = sample.src.x + 10;
    sample.dst.y = sample.src.y + 10;
    let sample_bv = get_bounding_volume(&sample);

    let nearest_a = tree.nearest(
        &sample,
        cmp_points,
        get_bounding_volume,
        Cutter { cut_limit: 1000, },
        bv_to_cut_point_sq_dist,
        bv_to_bv_sq_dist,
    ).next().unwrap().unwrap();

    let nearest_b = tree.iter()
        .flat_map(move |node| node.shapes())
        .map(|(shape, fragment)| (bv_to_bv_sq_dist(&sample_bv, fragment), shape, fragment))
        .min_by_key(|t| t.0)
        .unwrap();

    assert_eq!(nearest_a.dist, nearest_b.0);
    assert_eq!(nearest_a.shape, nearest_b.1);
    assert_eq!(nearest_a.shape_fragment, nearest_b.2);
}

#[test]
fn kdv_tree_nearest_big_o() {
    fn avg_dist_count(shapes_count: usize, avg_count: usize) -> (usize, usize, usize) {
        let tree = KdvTree::build(
            vec![Axis::X, Axis::Y],
            (0 .. shapes_count).map(|_| random_line(&mut rand::thread_rng())),
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            Cutter { cut_limit: 1000, }
        ).unwrap();
        let mut total = 0;
        for _ in 0 .. avg_count {
            let mut dist_cp_counter = 0;
            let mut dist_bv_counter = 0;
            tree.nearest(
                &random_line(&mut rand::thread_rng()),
                cmp_points,
                get_bounding_volume,
                Cutter { cut_limit: 1000, },
                |axis: &_, bounding_volume: &_, cut_point: &_| {
                    dist_cp_counter += 1;
                    bv_to_cut_point_sq_dist(axis, bounding_volume, cut_point)
                },
                |bv_a: &_, bv_b: &_| {
                    dist_bv_counter += 1;
                    bv_to_bv_sq_dist(bv_a, bv_b)
                },
            ).next().unwrap().unwrap();
            total += dist_cp_counter + dist_bv_counter;
        }
        let avg_total = total / avg_count;
        let mut max_depth = 0;
        let mut nodes_count = 0;
        for node in tree.iter() {
            if node.depth() > max_depth { max_depth = node.depth() }
            nodes_count += 1;
        }
        (avg_total, max_depth, nodes_count)
    }

    let stats_100 = avg_dist_count(100, 5);
    assert!(stats_100.0 * 10 < stats_100.2);
    let stats_1k = avg_dist_count(1000, 5);
    assert!(stats_1k.0 * 1000 < stats_1k.2);
    assert!(stats_1k.0 < stats_100.0 * 10);
    let stats_10k = avg_dist_count(10000, 5);
    assert!(stats_10k.0 * 10000 < stats_10k.2);
    assert!(stats_10k.0 < stats_1k.0 * 10);
}
