use std::cmp::{min, max, Ordering};
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

fn get_cut_point(_cut_axis: &Axis, points: &mut Iterator<Item = Point2d>) -> Option<Point2d> {
    let mut total = 0;
    let mut point_sum = Point2d { x: 0, y: 0, };
    for p in points {
        point_sum.x += p.x;
        point_sum.y += p.y;
        total += 1;
    }
    if total == 0 {
        None
    } else {
        Some(Point2d {
            x: (point_sum.x as f64 / total as f64) as i32,
            y: (point_sum.y as f64 / total as f64) as i32,
        })
    }
}

fn cut_volume(shape: &Line2d, fragment: &Rect2d, cut_axis: &Axis, cut_point: &Point2d) ->
    Result<Option<(Rect2d, Rect2d)>, ()>
{
    match cut_axis {
        &Axis::X => if cut_point.x >= fragment.lt.x && cut_point.x <= fragment.rb.x {
            if fragment.rb.x - fragment.lt.x < 10 {
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
            if fragment.rb.y - fragment.lt.y < 10 {
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

fn bv_to_cut_point_distance(axis: &Axis, bounding_volume: &Rect2d, cut_point: &Point2d) -> i64 {
    match axis {
        &Axis::X =>
            min((bounding_volume.lt.x - cut_point.x).abs() as i64, (bounding_volume.rb.x - cut_point.x).abs() as i64),
        &Axis::Y =>
            min((bounding_volume.lt.y - cut_point.y).abs() as i64, (bounding_volume.rb.y - cut_point.y).abs() as i64),
    }
}

fn bv_to_bv_distance(bounding_volume_a: &Rect2d, bounding_volume_b: &Rect2d) -> i64 {
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

#[test]
fn kdv_tree_basic() {
    let shapes = vec![Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, }];
    let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, cmp_points, get_bounding_volume, get_cut_point, cut_volume).unwrap();

    assert_eq!(
        tree.intersects(
            &Line2d { src: Point2d { x: 116, y: 116, }, dst: Point2d { x: 180, y: 180, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            cut_volume,
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
            cut_volume,
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
            cut_volume,
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
            cut_volume,
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
            cut_volume,
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
    let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, cmp_points, get_bounding_volume, get_cut_point, cut_volume).unwrap();

    assert_eq!(
        tree.intersects(
            &Line2d { src: Point2d { x: 70, y: 45, }, dst: Point2d { x: 75, y: 50, }, },
            cmp_points,
            get_bounding_volume,
            get_cut_point,
            cut_volume,
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
            cut_volume,
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
            cut_volume,
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
    let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, cmp_points, get_bounding_volume, get_cut_point, cut_volume).unwrap();

    let nearest: Vec<_> = tree
        .nearest(
            &Line2d { src: Point2d { x: 44, y: 8, }, dst: Point2d { x: 52, y: 12, }, },
            cmp_points,
            get_bounding_volume,
            bv_to_cut_point_distance,
            bv_to_bv_distance,
        )
        .collect();
    assert_eq!(nearest, vec![
        NearestShape { dist: 61, shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } } },
        NearestShape { dist: 724, shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } } },
        NearestShape { dist: 1553, shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } } },
    ]);
    let nearest: Vec<_> = tree
        .nearest(
            &Line2d { src: Point2d { x: 44, y: 88, }, dst: Point2d { x: 52, y: 92, }, },
            cmp_points,
            get_bounding_volume,
            bv_to_cut_point_distance,
            bv_to_bv_distance,
        )
        .collect();
    assert_eq!(nearest, vec![
        NearestShape { dist: 884, shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } } },
        NearestShape { dist: 2180, shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } } },
        NearestShape { dist: 5501, shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } } },
    ]);
    let nearest: Vec<_> = tree
        .nearest(
            &Line2d { src: Point2d { x: 144, y: 48, }, dst: Point2d { x: 152, y: 52, }, },
            cmp_points,
            get_bounding_volume,
            bv_to_cut_point_distance,
            bv_to_bv_distance,
        )
        .collect();
    assert_eq!(nearest, vec![
        NearestShape { dist: 4660, shape: &Line2d { src: Point2d { x: 80, y: 16 }, dst: Point2d { x: 80, y: 80 } } },
        NearestShape { dist: 5770, shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 80 } } },
        NearestShape { dist: 6197, shape: &Line2d { src: Point2d { x: 16, y: 16 }, dst: Point2d { x: 80, y: 16 } } },
    ]);
}
