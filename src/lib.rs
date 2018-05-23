use std::iter;
use std::cmp::Ordering;

pub trait Axis<Point> {
    fn cmp_points(&self, a: &Point, b: &Point) -> Ordering;
}

pub trait BoundingVolume {
    type Point;

    fn min_corner(&self) -> Self::Point;
    fn max_corner(&self) -> Self::Point;
}

pub trait VolumeManager<S, A> {
    type BoundingVolume: BoundingVolume;
    type Error;

    fn bounding_volume(&self, shape: &S) -> Self::BoundingVolume;

    fn cut_point<I>(&mut self, cut_axis: &A, points: I) ->
        Option<<Self::BoundingVolume as BoundingVolume>::Point>
        where I: Iterator<Item = <Self::BoundingVolume as BoundingVolume>::Point>;

    fn cut(
        &mut self,
        shape: &S,
        fragment: &Self::BoundingVolume,
        cut_axis: &A,
        cut_point: &<Self::BoundingVolume as BoundingVolume>::Point,
    ) -> Result<Option<(Self::BoundingVolume, Self::BoundingVolume)>, Self::Error>;
}

pub struct KdvTree<A, P, B, S> {
    axis: Vec<A>,
    shapes: Vec<S>,
    root: KdvNode<P, B>,
}

impl<A, P, B, S> KdvTree<A, P, B, S>
    where A: Axis<P>,
          B: BoundingVolume<Point = P>,
{
    pub fn build<IA, II, M>(axis_it: IA, shapes_it: II, manager: &mut M) -> Result<KdvTree<A, P, B, S>, M::Error>
        where IA: IntoIterator<Item = A>,
              II: IntoIterator<Item = S>,
              M: VolumeManager<S, A, BoundingVolume = B>,
    {
        let axis: Vec<_> = axis_it.into_iter().collect();
        let shapes: Vec<_> = shapes_it.into_iter().collect();
        let root_shapes: Vec<_> = shapes
            .iter()
            .enumerate()
            .map(|(i, s)| ShapeFragment {
                bounding_volume: manager.bounding_volume(s),
                shape_id: i,
            })
            .collect();
        Ok(KdvTree {
            root: KdvNode::build(0, &axis, &shapes, root_shapes, manager)?,
            axis, shapes,
        })
    }

    pub fn intersects<'t, 's, 'm, SN, M>(&'t self, shape: &'s SN, manager: &'m mut M) ->
        IntersectIter<'t, 's, 'm, A, P, S, B, SN, M::BoundingVolume, M>
        where M: VolumeManager<SN, A>,
    {
        IntersectIter {
            needle: shape,
            axis: &self.axis,
            shapes: &self.shapes,
            queue: vec![TraverseTask::Explore {
                node: &self.root,
                needle_fragment: manager.bounding_volume(shape),
            }],
            manager,
        }
    }
}

struct ShapeFragment<B> {
    bounding_volume: B,
    shape_id: usize,
}

struct KdvNode<P, B> {
    shapes: Vec<ShapeFragment<B>>,
    children: KdvNodeChildren<P, B>,
}

struct CutPoint<P> {
    axis_index: usize,
    point: P,
}

enum KdvNodeChildren<P, B> {
    Missing,
    OnlyLeft { cut_point: CutPoint<P>, child: Box<KdvNode<P, B>>, },
    OnlyRight { cut_point: CutPoint<P>, child: Box<KdvNode<P, B>>, },
    Both { cut_point: CutPoint<P>, left: Box<KdvNode<P, B>>, right: Box<KdvNode<P, B>>, },
}

impl<P, B> KdvNode<P, B> {
    fn build<A, S, M>(depth: usize, axis: &[A], shapes: &[S], mut node_shapes: Vec<ShapeFragment<B>>, manager: &mut M) ->
        Result<KdvNode<P, B>, M::Error>
        where B: BoundingVolume<Point = P>,
              A: Axis<P>,
              M: VolumeManager<S, A, BoundingVolume = B>,
    {
        // locate cut point for coords
        let cut_axis_index = depth % axis.len();
        let cut_axis = &axis[cut_axis_index];
        let maybe_cut_point = manager.cut_point(
            cut_axis,
            node_shapes
                .iter()
                .flat_map(|sf| {
                    let bvol = &sf.bounding_volume;
                    iter::once(bvol.min_corner())
                        .chain(iter::once(bvol.max_corner()))
                }),
        );

        if let Some(cut_point) = maybe_cut_point {
            // distribute shapes among children
            let mut left_shapes = Vec::new();
            let mut right_shapes = Vec::new();
            let mut head = 0;
            while head < node_shapes.len() {
                let ShapeFragment { shape_id, bounding_volume, } = node_shapes.swap_remove(head);
                let owner = shape_owner(&shapes[shape_id], bounding_volume, cut_axis, &cut_point, manager)?;
                match owner {
                    ShapeOwner::Me(bounding_volume) => {
                        let tail = node_shapes.len();
                        node_shapes.push(ShapeFragment { shape_id, bounding_volume, });
                        node_shapes.swap(head, tail);
                        head += 1;
                    },
                    ShapeOwner::Left(bounding_volume) =>
                        left_shapes.push(ShapeFragment { shape_id, bounding_volume, }),
                    ShapeOwner::Right(bounding_volume) =>
                        right_shapes.push(ShapeFragment { shape_id, bounding_volume, }),
                    ShapeOwner::Both { left_bvol, right_bvol, } => {
                        left_shapes.push(ShapeFragment { shape_id, bounding_volume: left_bvol, });
                        right_shapes.push(ShapeFragment { shape_id, bounding_volume: right_bvol, });
                    },
                }
            }
            // construct the node
            let children = if left_shapes.is_empty() && right_shapes.is_empty() {
                KdvNodeChildren::Missing
            } else if left_shapes.is_empty() {
                KdvNodeChildren::OnlyRight {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes, manager)?),
                }
            } else if right_shapes.is_empty() {
                KdvNodeChildren::OnlyLeft {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes, manager)?),
                }
            } else {
                KdvNodeChildren::Both {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    left: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes, manager)?),
                    right: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes, manager)?),
                }
            };
            Ok(KdvNode { shapes: node_shapes, children, })
        } else {
            // no cut point choosen, keep all shapes in current node
            Ok(KdvNode {
                shapes: node_shapes,
                children: KdvNodeChildren::Missing,
            })
        }
    }
}

enum ShapeOwner<B> {
    Me(B),
    Left(B),
    Right(B),
    Both { left_bvol: B, right_bvol: B, },
}

fn shape_owner<A, P, B, S, M>(shape: &S, fragment: B, cut_axis: &A, cut_point: &P, manager: &mut M) ->
    Result<ShapeOwner<B>, M::Error>
    where A: Axis<P>,
          B: BoundingVolume<Point = P>,
          M: VolumeManager<S, A, BoundingVolume = B>,
{
    let min_corner = fragment.min_corner();
    let max_corner = fragment.max_corner();
    Ok(match (cut_axis.cmp_points(&min_corner, cut_point), cut_axis.cmp_points(&max_corner, cut_point)) {
        (Ordering::Less, Ordering::Less) | (Ordering::Less, Ordering::Equal) =>
            ShapeOwner::Left(fragment),
        (Ordering::Greater, Ordering::Greater) | (Ordering::Equal, Ordering::Greater) =>
            ShapeOwner::Right(fragment),
        _ => if let Some((left_bvol, right_bvol)) = manager.cut(shape, &fragment, cut_axis, cut_point)? {
            ShapeOwner::Both { left_bvol, right_bvol, }
        } else {
            ShapeOwner::Me(fragment)
        }
    })
}

enum TraverseTask<'t, P: 't, BS: 't, BN> {
    Explore { node: &'t KdvNode<P, BS>, needle_fragment: BN, },
    Intersect { needle_fragment: BN, shape_fragment: &'t ShapeFragment<BS>, axis_counter: usize, },
}

pub struct IntersectIter<'t, 's, 'm, A: 't, P: 't, SS: 't, BS: 't, SN: 's, BN, M: 'm> {
    needle: &'s SN,
    axis: &'t [A],
    shapes: &'t [SS],
    queue: Vec<TraverseTask<'t, P, BS, BN>>,
    manager: &'m mut M,
}

#[derive(Clone, PartialEq, Debug)]
pub struct Intersection<'t, SS: 't, BS: 't, BN> {
    pub shape: &'t SS,
    pub shape_fragment: &'t BS,
    pub needle_fragment: BN,
}

impl<'t, 's, 'm, A, P, SS, BS, SN, BN, M> Iterator for IntersectIter<'t, 's, 'm, A, P, SS, BS, SN, BN, M>
    where A: Axis<P>,
          BS: BoundingVolume<Point = P>,
          BN: BoundingVolume<Point = BS::Point> + Clone,
          M: VolumeManager<SN, A, BoundingVolume = BN>,
{
    type Item = Result<Intersection<'t, SS, BS, BN>, M::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        'outer: while let Some(task) = self.queue.pop() {
            match task {
                TraverseTask::Explore { node, needle_fragment, } => {
                    match node.children {
                        KdvNodeChildren::Missing =>
                            (),
                        KdvNodeChildren::OnlyLeft { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, self.manager) {
                                Ok(ShapeOwner::Me(needle_fragment)) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                Ok(ShapeOwner::Left(needle_fragment)) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                Ok(ShapeOwner::Right(..)) =>
                                    (),
                                Ok(ShapeOwner::Both { left_bvol: needle_fragment, .. }) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                Err(error) => {
                                    self.queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        },
                        KdvNodeChildren::OnlyRight { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, self.manager) {
                                Ok(ShapeOwner::Me(needle_fragment)) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                Ok(ShapeOwner::Left(..)) =>
                                    (),
                                Ok(ShapeOwner::Right(needle_fragment)) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                Ok(ShapeOwner::Both { right_bvol: needle_fragment, .. }) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                Err(error) => {
                                    self.queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        },
                        KdvNodeChildren::Both { ref cut_point, ref left, ref right, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, self.manager) {
                                Ok(ShapeOwner::Me(fragment)) => {
                                    self.queue.push(TraverseTask::Explore { node: left, needle_fragment: fragment.clone(), });
                                    self.queue.push(TraverseTask::Explore { node: right, needle_fragment: fragment, });
                                },
                                Ok(ShapeOwner::Left(needle_fragment)) =>
                                    self.queue.push(TraverseTask::Explore { node: left, needle_fragment, }),
                                Ok(ShapeOwner::Right(needle_fragment)) =>
                                    self.queue.push(TraverseTask::Explore { node: right, needle_fragment, }),
                                Ok(ShapeOwner::Both { left_bvol, right_bvol, }) => {
                                    self.queue.push(TraverseTask::Explore { node: left, needle_fragment: left_bvol, });
                                    self.queue.push(TraverseTask::Explore { node: right, needle_fragment: right_bvol, });
                                },
                                Err(error) => {
                                    self.queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        },
                    }
                    for shape_fragment in node.shapes.iter() {
                        self.queue.push(TraverseTask::Intersect {
                            shape_fragment,
                            needle_fragment: needle_fragment.clone(),
                            axis_counter: 0,
                        });
                    }
                },
                TraverseTask::Intersect { shape_fragment, needle_fragment, mut axis_counter, } => {
                    let no_intersection = self.axis.iter().any(|axis| {
                        let needle_min = needle_fragment.min_corner();
                        let needle_max = needle_fragment.max_corner();
                        let shape_min = shape_fragment.bounding_volume.min_corner();
                        let shape_max = shape_fragment.bounding_volume.max_corner();
                        (axis.cmp_points(&needle_min, &shape_max) == Ordering::Greater ||
                         axis.cmp_points(&needle_max, &shape_min) == Ordering::Less)
                    });
                    if no_intersection {
                        continue;
                    }
                    let axis_total = self.axis.len();
                    for _ in 0 .. axis_total {
                        let cut_axis = &self.axis[axis_counter % axis_total];
                        axis_counter += 1;
                        let maybe_cut_point = self.manager.cut_point(
                            cut_axis,
                            iter::once(needle_fragment.min_corner())
                                .chain(iter::once(needle_fragment.max_corner())),
                        );
                        if let Some(cut_point) = maybe_cut_point {
                            match self.manager.cut(self.needle, &needle_fragment, cut_axis, &cut_point) {
                                Ok(Some((left_fragment, right_fragment))) => {
                                    self.queue.push(TraverseTask::Intersect {
                                        shape_fragment: shape_fragment.clone(),
                                        needle_fragment: left_fragment,
                                        axis_counter,
                                    });
                                    self.queue.push(TraverseTask::Intersect {
                                        shape_fragment: shape_fragment,
                                        needle_fragment: right_fragment,
                                        axis_counter,
                                    });
                                    continue 'outer;
                                },
                                Ok(None) =>
                                    (),
                                Err(error) => {
                                    self.queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        }
                    }
                    return Some(Ok(Intersection {
                        shape: &self.shapes[shape_fragment.shape_id],
                        shape_fragment: &shape_fragment.bounding_volume,
                        needle_fragment,
                    }));
                },
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::{min, max, Ordering};
    use super::{KdvTree, Intersection};

    #[derive(Clone, Copy, PartialEq, Debug)]
    struct Point2d { x: i32, y: i32, }

    #[derive(Clone, Debug)]
    enum Axis { X, Y, }

    impl super::Axis<Point2d> for Axis {
        fn cmp_points(&self, a: &Point2d, b: &Point2d) -> Ordering {
            match self {
                &Axis::X => a.x.cmp(&b.x),
                &Axis::Y => a.y.cmp(&b.y),
            }
        }
    }

    #[derive(PartialEq, Clone, Debug)]
    struct Rect2d { lt: Point2d, rb: Point2d, }

    impl super::BoundingVolume for Rect2d {
        type Point = Point2d;

        fn min_corner(&self) -> Self::Point { self.lt }
        fn max_corner(&self) -> Self::Point { self.rb }
    }

    #[derive(PartialEq, Debug)]
    struct Line2d { src: Point2d, dst: Point2d, }

    struct Manager;

    impl super::VolumeManager<Line2d, Axis> for Manager {
        type BoundingVolume = Rect2d;
        type Error = ();

        fn bounding_volume(&self, shape: &Line2d) -> Self::BoundingVolume {
            Rect2d {
                lt: Point2d { x: min(shape.src.x, shape.dst.x), y: min(shape.src.y, shape.dst.y), },
                rb: Point2d { x: max(shape.src.x, shape.dst.x), y: max(shape.src.y, shape.dst.y), },
            }
        }

        fn cut_point<I>(&mut self, _cut_axis: &Axis, points: I) -> Option<Point2d> where I: Iterator<Item = Point2d> {
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

        fn cut(&mut self, shape: &Line2d, fragment: &Rect2d, cut_axis: &Axis, cut_point: &Point2d) ->
            Result<Option<(Rect2d, Rect2d)>, Self::Error>
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
    }

    #[test]
    fn kdv_tree_basic() {
        let shapes = vec![Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, }];
        let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, &mut Manager).unwrap();

        assert_eq!(
            tree.intersects(&Line2d { src: Point2d { x: 116, y: 116, }, dst: Point2d { x: 180, y: 180, }, }, &mut Manager)
                .collect::<Result<Vec<_>, _>>(),
            Ok(vec![])
        );
        assert_eq!(
            tree.intersects(&Line2d { src: Point2d { x: 32, y: 48, }, dst: Point2d { x: 48, y: 64, }, }, &mut Manager)
                .collect::<Result<Vec<_>, _>>(),
            Ok(vec![])
        );
        assert_eq!(
            tree.intersects(&Line2d { src: Point2d { x: 48, y: 32, }, dst: Point2d { x: 64, y: 48, }, }, &mut Manager)
                .collect::<Result<Vec<_>, _>>(),
            Ok(vec![])
        );

        let intersects: Result<Vec<_>, _> = tree
            .intersects(&Line2d { src: Point2d { x: 16, y: 64, }, dst: Point2d { x: 80, y: 64, }, }, &mut Manager)
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
            .intersects(&Line2d { src: Point2d { x: 64, y: 16, }, dst: Point2d { x: 64, y: 80, }, }, &mut Manager)
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
        let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes, &mut Manager).unwrap();

        assert_eq!(
            tree.intersects(&Line2d { src: Point2d { x: 70, y: 45, }, dst: Point2d { x: 75, y: 50, }, }, &mut Manager)
                .collect::<Result<Vec<_>, _>>(),
            Ok(vec![])
        );

        let intersects: Result<Vec<_>, _> = tree
            .intersects(&Line2d { src: Point2d { x: 8, y: 48, }, dst: Point2d { x: 88, y: 48, }, }, &mut Manager)
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
            .intersects(&Line2d { src: Point2d { x: 40, y: 10, }, dst: Point2d { x: 90, y: 60, }, }, &mut Manager)
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
}
