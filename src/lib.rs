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

pub trait Shape {
    type BoundingBox: BoundingBox;

    fn bounding_box(&self) -> Self::BoundingBox;
    fn cut<A>(
        &self,
        fragment: &Self::BoundingBox,
        cut_axis: &A,
        cut_point: &<Self::BoundingBox as BoundingBox>::Point
    ) -> Option<(Self::BoundingBox, Self::BoundingBox)>
        where A: Axis<<Self::BoundingBox as BoundingBox>::Point>;
}

pub struct KdvTree<A, P, B, S> {
    axis: Vec<A>,
    shapes: Vec<S>,
    root: KdvNode<P, B>,
}

impl<A, P, B, S> KdvTree<A, P, B, S>
    where A: Axis<P>,
          B: BoundingBox<Point = P>,
          S: Shape<BoundingBox = B>,
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
        where S: Shape<BoundingBox = B>,
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
          S: Shape<BoundingBox = B>,
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

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
