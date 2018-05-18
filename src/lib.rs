use std::iter;
use std::cmp::Ordering;

pub trait Axis<Point> {
    fn cut_point<I>(&self, points: I) -> Option<Point> where I: Iterator<Item = Point>;
    fn cmp_points(&self, a: &Point, b: &Point) -> Ordering;
}

pub trait BoundingBox {
    type Point;

    fn min_corner(&self) -> Self::Point;
    fn max_corner(&self) -> Self::Point;
}

pub trait Shape<A> {
    type BoundingBox: BoundingBox;

    fn bounding_box(&self) -> Self::BoundingBox;
    fn cut(
        &self,
        fragment: &Self::BoundingBox,
        cut_axis: &A,
        cut_point: &<Self::BoundingBox as BoundingBox>::Point
    ) -> Option<(Self::BoundingBox, Self::BoundingBox)>;
}

pub struct KdvTree<A, P, B, S> {
    axis: Vec<A>,
    shapes: Vec<S>,
    root: KdvNode<P, B>,
}

impl<A, P, B, S> KdvTree<A, P, B, S>
    where A: Axis<P>,
          B: BoundingBox<Point = P>,
          S: Shape<A, BoundingBox = B>,
{
    pub fn build<IA, II>(axis_it: IA, shapes_it: II) -> KdvTree<A, P, B, S>
        where IA: IntoIterator<Item = A>,
              II: IntoIterator<Item = S>,
    {
        let axis: Vec<_> = axis_it.into_iter().collect();
        let shapes: Vec<_> = shapes_it.into_iter().collect();
        let root_shapes: Vec<_> = shapes
            .iter()
            .enumerate()
            .map(|(i, s)| ShapeFragment {
                bounding_box: s.bounding_box(),
                shape_id: i,
            })
            .collect();
        KdvTree {
            root: KdvNode::build(0, &axis, &shapes, root_shapes),
            axis, shapes,
        }
    }

    pub fn intersects<'t, 's, SN>(&'t self, shape: &'s SN) -> IntersectIter<'t, 's, A, P, S, B, SN, SN::BoundingBox>
        where SN: Shape<A, BoundingBox = S::BoundingBox>
    {
        IntersectIter {
            needle: shape,
            axis: &self.axis,
            shapes: &self.shapes,
            queue: vec![TraverseTask::Explore {
                node: &self.root,
                needle_fragment: shape.bounding_box(),
            }],
        }
    }
}

struct ShapeFragment<B> {
    bounding_box: B,
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
    fn build<A, S>(depth: usize, axis: &[A], shapes: &[S], mut node_shapes: Vec<ShapeFragment<B>>) -> KdvNode<P, B>
        where S: Shape<A, BoundingBox = B>,
              B: BoundingBox<Point = P>,
              A: Axis<P>,
    {
        // locate cut point for coords
        let cut_axis_index = depth % axis.len();
        let cut_axis = &axis[cut_axis_index];
        let maybe_cut_point = cut_axis.cut_point(
            node_shapes
                .iter()
                .flat_map(|sf| {
                    let bbox = &sf.bounding_box;
                    iter::once(bbox.min_corner())
                        .chain(iter::once(bbox.max_corner()))
                })
        );

        if let Some(cut_point) = maybe_cut_point {
            // distribute shapes among children
            let mut left_shapes = Vec::new();
            let mut right_shapes = Vec::new();
            let mut head = 0;
            while head < node_shapes.len() {
                let ShapeFragment { shape_id, bounding_box, } = node_shapes.swap_remove(head);
                let owner = shape_owner(&shapes[shape_id], bounding_box, cut_axis, &cut_point);
                match owner {
                    ShapeOwner::Me(bounding_box) => {
                        let tail = node_shapes.len();
                        node_shapes.push(ShapeFragment { shape_id, bounding_box, });
                        node_shapes.swap(head, tail);
                        head += 1;
                    },
                    ShapeOwner::Left(bounding_box) =>
                        left_shapes.push(ShapeFragment { shape_id, bounding_box, }),
                    ShapeOwner::Right(bounding_box) =>
                        right_shapes.push(ShapeFragment { shape_id, bounding_box, }),
                    ShapeOwner::Both { left_bbox, right_bbox, } => {
                        left_shapes.push(ShapeFragment { shape_id, bounding_box: left_bbox, });
                        right_shapes.push(ShapeFragment { shape_id, bounding_box: right_bbox, });
                    },
                }
            }
            // construct the node
            let children = if left_shapes.is_empty() && right_shapes.is_empty() {
                KdvNodeChildren::Missing
            } else if left_shapes.is_empty() {
                KdvNodeChildren::OnlyRight {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes)),
                }
            } else if right_shapes.is_empty() {
                KdvNodeChildren::OnlyLeft {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes)),
                }
            } else {
                KdvNodeChildren::Both {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    left: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes)),
                    right: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes)),
                }
            };
            KdvNode { shapes: node_shapes, children, }
        } else {
            // no cut point choosen, keep all shapes in current node
            KdvNode {
                shapes: node_shapes,
                children: KdvNodeChildren::Missing,
            }
        }
    }
}

enum ShapeOwner<B> {
    Me(B),
    Left(B),
    Right(B),
    Both { left_bbox: B, right_bbox: B, },
}

fn shape_owner<A, P, B, S>(shape: &S, fragment: B, cut_axis: &A, cut_point: &P) -> ShapeOwner<B>
    where A: Axis<P>,
          B: BoundingBox<Point = P>,
          S: Shape<A, BoundingBox = B>,
{
    let min_corner = fragment.min_corner();
    let max_corner = fragment.max_corner();
    match (cut_axis.cmp_points(&min_corner, cut_point), cut_axis.cmp_points(&max_corner, cut_point)) {
        (Ordering::Less, Ordering::Less) | (Ordering::Less, Ordering::Equal) =>
            ShapeOwner::Left(fragment),
        (Ordering::Greater, Ordering::Greater) | (Ordering::Equal, Ordering::Greater) =>
            ShapeOwner::Right(fragment),
        _ => if let Some((left_bbox, right_bbox)) = shape.cut(&fragment, cut_axis, cut_point) {
            ShapeOwner::Both { left_bbox, right_bbox, }
        } else {
            ShapeOwner::Me(fragment)
        }
    }
}

enum TraverseTask<'t, P: 't, BS: 't, BN> {
    Explore { node: &'t KdvNode<P, BS>, needle_fragment: BN, },
    Intersect { needle_fragment: BN, shape_fragment: &'t ShapeFragment<BS>, axis_counter: usize, },
}

pub struct IntersectIter<'t, 's, A: 't, P: 't, SS: 't, BS: 't, SN: 's, BN> {
    needle: &'s SN,
    axis: &'t [A],
    shapes: &'t [SS],
    queue: Vec<TraverseTask<'t, P, BS, BN>>,
}

#[derive(Clone, PartialEq, Debug)]
pub struct Intersection<'t, SS: 't, BS: 't, BN> {
    pub shape: &'t SS,
    pub shape_fragment: &'t BS,
    pub needle_fragment: BN,
}

impl<'t, 's, A, P, SS, BS, SN, BN> Iterator for IntersectIter<'t, 's, A, P, SS, BS, SN, BN>
    where A: Axis<P>,
          SS: Shape<A, BoundingBox = BS>,
          BS: BoundingBox<Point = P>,
          SN: Shape<A, BoundingBox = BN>,
          BN: BoundingBox<Point = BS::Point> + Clone,
{
    type Item = Intersection<'t, SS, BS, BN>;

    fn next(&mut self) -> Option<Self::Item> {
        'outer: while let Some(task) = self.queue.pop() {
            match task {
                TraverseTask::Explore { node, needle_fragment, } => {
                    match node.children {
                        KdvNodeChildren::Missing =>
                            (),
                        KdvNodeChildren::OnlyLeft { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point) {
                                ShapeOwner::Me(needle_fragment) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                ShapeOwner::Left(needle_fragment) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                ShapeOwner::Right(..) =>
                                    (),
                                ShapeOwner::Both { left_bbox: needle_fragment, .. } =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                            }
                        },
                        KdvNodeChildren::OnlyRight { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point) {
                                ShapeOwner::Me(needle_fragment) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                ShapeOwner::Left(..) =>
                                    (),
                                ShapeOwner::Right(needle_fragment) =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                                ShapeOwner::Both { right_bbox: needle_fragment, .. } =>
                                    self.queue.push(TraverseTask::Explore { node: child, needle_fragment, }),
                            }
                        },
                        KdvNodeChildren::Both { ref cut_point, ref left, ref right, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point) {
                                ShapeOwner::Me(fragment) => {
                                    self.queue.push(TraverseTask::Explore { node: left, needle_fragment: fragment.clone(), });
                                    self.queue.push(TraverseTask::Explore { node: right, needle_fragment: fragment, });
                                },
                                ShapeOwner::Left(needle_fragment) =>
                                    self.queue.push(TraverseTask::Explore { node: left, needle_fragment, }),
                                ShapeOwner::Right(needle_fragment) =>
                                    self.queue.push(TraverseTask::Explore { node: right, needle_fragment, }),
                                ShapeOwner::Both { left_bbox, right_bbox, } => {
                                    self.queue.push(TraverseTask::Explore { node: left, needle_fragment: left_bbox, });
                                    self.queue.push(TraverseTask::Explore { node: right, needle_fragment: right_bbox, });
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
                        let shape_min = shape_fragment.bounding_box.min_corner();
                        let shape_max = shape_fragment.bounding_box.max_corner();
                        (axis.cmp_points(&needle_min, &shape_max) == Ordering::Greater ||
                         axis.cmp_points(&needle_max, &shape_min) == Ordering::Less)
                    });
                    if no_intersection {
                        continue;
                    }
                    let axis_total = self.axis.len();
                    for _ in 0 .. axis_total {
                        let cut_axis = &self.axis[axis_counter % axis_total];
                        let maybe_cut_point = cut_axis.cut_point(
                            iter::once(needle_fragment.min_corner())
                                .chain(iter::once(needle_fragment.max_corner()))
                        );
                        if let Some(cut_point) = maybe_cut_point {
                            if let Some((left_fragment, right_fragment)) = self.needle.cut(&needle_fragment, cut_axis, &cut_point) {
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
                            }
                        }
                        axis_counter += 1;
                    }
                    return Some(Intersection {
                        shape: &self.shapes[shape_fragment.shape_id],
                        shape_fragment: &shape_fragment.bounding_box,
                        needle_fragment,
                    });
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
        fn cut_point<I>(&self, points: I) -> Option<Point2d> where I: Iterator<Item = Point2d> {
            let mut total = 0;
            let mut sum = 0;
            for p in points {
                sum += match self {
                    &Axis::X => p.x,
                    &Axis::Y => p.y,
                };
                total += 1;
            }
            if total == 0 {
                None
            } else {
                let mid = sum / total;
                Some(match self {
                    &Axis::X => Point2d { x: mid, y: 0, },
                    &Axis::Y => Point2d { x: 0, y: mid, },
                })
            }
        }

        fn cmp_points(&self, a: &Point2d, b: &Point2d) -> Ordering {
            match self {
                &Axis::X => a.x.cmp(&b.x),
                &Axis::Y => a.y.cmp(&b.y),
            }
        }
    }

    #[derive(PartialEq, Clone, Debug)]
    struct Rect2d { lt: Point2d, rb: Point2d, }

    impl super::BoundingBox for Rect2d {
        type Point = Point2d;

        fn min_corner(&self) -> Self::Point { self.lt }
        fn max_corner(&self) -> Self::Point { self.rb }
    }

    #[derive(PartialEq, Debug)]
    struct Line2d { src: Point2d, dst: Point2d, }

    impl super::Shape<Axis> for Line2d {
        type BoundingBox = Rect2d;

        fn bounding_box(&self) -> Self::BoundingBox {
            Rect2d {
                lt: Point2d { x: min(self.src.x, self.dst.x), y: min(self.src.y, self.dst.y), },
                rb: Point2d { x: max(self.src.x, self.dst.x), y: max(self.src.y, self.dst.y), },
            }
        }

        fn cut(
            &self,
            fragment: &Self::BoundingBox,
            cut_axis: &Axis,
            cut_point: &<Self::BoundingBox as super::BoundingBox>::Point
        )
            -> Option<(Self::BoundingBox, Self::BoundingBox)>
        {
            let bbox = self.bounding_box();
            let (side, x, y) = match cut_axis {
                &Axis::X => if cut_point.x >= fragment.lt.x && cut_point.x <= fragment.rb.x {
                    let factor = (cut_point.x - bbox.lt.x) as f64 / (bbox.rb.x - bbox.lt.x) as f64;
                    (fragment.rb.x - fragment.lt.x, cut_point.x, bbox.lt.y + (factor * (bbox.rb.y - bbox.lt.y) as f64) as i32)
                } else {
                    return None;
                },
                &Axis::Y => if cut_point.y >= fragment.lt.y && cut_point.y <= fragment.rb.y {
                    let factor = (cut_point.y - bbox.lt.y) as f64 / (bbox.rb.y - bbox.lt.y) as f64;
                    (fragment.rb.y - fragment.lt.y, bbox.lt.x + (factor * (bbox.rb.x - bbox.lt.x) as f64) as i32, cut_point.y)
                } else {
                    return None;
                },
            };
            if side < 10 {
                None
            } else {
                Some((Rect2d { lt: fragment.lt, rb: Point2d { x, y, } }, Rect2d { lt: Point2d { x, y, }, rb: fragment.rb, }))
            }
        }
    }

    #[test]
    fn kdv_tree_basic() {
        let shapes = vec![Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, }];
        let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes);

        assert_eq!(tree.intersects(&Line2d { src: Point2d { x: 116, y: 116, }, dst: Point2d { x: 180, y: 180, }, }).collect::<Vec<_>>(), vec![]);
        assert_eq!(tree.intersects(&Line2d { src: Point2d { x: 32, y: 48, }, dst: Point2d { x: 48, y: 64, }, }).collect::<Vec<_>>(), vec![]);
        assert_eq!(tree.intersects(&Line2d { src: Point2d { x: 48, y: 32, }, dst: Point2d { x: 64, y: 48, }, }).collect::<Vec<_>>(), vec![]);

        let intersects: Vec<_> = tree
            .intersects(&Line2d { src: Point2d { x: 16, y: 64, }, dst: Point2d { x: 80, y: 64, }, })
            .collect();
        assert_eq!(intersects, vec![
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
        ]);
        let intersects: Vec<_> = tree
            .intersects(&Line2d { src: Point2d { x: 64, y: 16, }, dst: Point2d { x: 64, y: 80, }, })
            .collect();
        assert_eq!(intersects, [
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
        ]);
    }

    #[test]
    fn kdv_tree_triangle() {
        let shapes = vec![
            Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 16, }, },
            Line2d { src: Point2d { x: 16, y: 16, }, dst: Point2d { x: 80, y: 80, }, },
            Line2d { src: Point2d { x: 80, y: 16, }, dst: Point2d { x: 80, y: 80, }, },
        ];
        let tree = KdvTree::build(vec![Axis::X, Axis::Y], shapes);

        assert_eq!(tree.intersects(&Line2d { src: Point2d { x: 70, y: 45, }, dst: Point2d { x: 75, y: 50, }, }).collect::<Vec<_>>(), vec![]);

        let intersects: Vec<_> = tree
            .intersects(&Line2d { src: Point2d { x: 8, y: 48, }, dst: Point2d { x: 88, y: 48, }, })
            .collect();
        assert_eq!(intersects, vec![
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
        ]);
        let intersects: Vec<_> = tree
            .intersects(&Line2d { src: Point2d { x: 40, y: 10, }, dst: Point2d { x: 90, y: 60, }, })
            .collect();
        assert_eq!(intersects, vec![
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
        ]);
    }
}
