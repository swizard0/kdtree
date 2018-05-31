use std::cmp::Ordering;
use std::iter::{self, FromIterator};
use std::collections::{BinaryHeap, HashSet};
use std::collections::binary_heap::PeekMut;

#[cfg(test)]
mod tests;

pub trait BoundingVolume<P> {
    fn min_corner(&self) -> P;
    fn max_corner(&self) -> P;
}

pub trait CmpPoints<A, P> {
    fn cmp_points(&self, axis: &A, a: &P, b: &P) -> Ordering;
}

impl<A, P, F> CmpPoints<A, P> for F where F: Fn(&A, &P, &P) -> Ordering {
    fn cmp_points(&self, axis: &A, a: &P, b: &P) -> Ordering {
        (self)(axis, a, b)
    }
}

pub trait GetBoundingVolume<B, S> {
    fn bounding_volume(&self, shape: &S) -> B;
}

impl<B, S, F> GetBoundingVolume<B, S> for F where F: Fn(&S) -> B {
    fn bounding_volume(&self, shape: &S) -> B {
        (self)(shape)
    }
}

pub trait GetCutPoint<A, P> {
    fn cut_point<I>(&mut self, cut_axis: &A, points: I) -> Option<P>
        where I: Iterator<Item = P>;
}

impl<A, P, F> GetCutPoint<A, P> for F where F: FnMut(&A, &mut Iterator<Item = P>) -> Option<P> {
    fn cut_point<I>(&mut self, cut_axis: &A, mut points: I) -> Option<P>
        where I: Iterator<Item = P>
    {
        (self)(cut_axis, &mut points)
    }
}

pub trait BoundingVolumesCutter<A, P, B, S> {
    type Error;

    fn cut(&mut self, shape: &S, fragment: &B, cut_axis: &A, cut_point: &P) -> Result<Option<(B, B)>, Self::Error>;
}

impl<A, P, B, S, F, E> BoundingVolumesCutter<A, P, B, S> for F where F: FnMut(&S, &B, &A, &P) -> Result<Option<(B, B)>, E> {
    type Error = E;

    fn cut(&mut self, shape: &S, fragment: &B, cut_axis: &A, cut_point: &P) -> Result<Option<(B, B)>, Self::Error> {
        (self)(shape, fragment, cut_axis, cut_point)
    }
}

pub trait DistanceBVCP<A, P, B> {
    type Dist;

    fn bv_to_cut_point_distance(&self, axis: &A, bounding_volume: &B, cut_point: &P) -> Self::Dist;
}

impl<A, P, B, F, D> DistanceBVCP<A, P, B> for F where F: Fn(&A, &B, &P) -> D {
    type Dist = D;

    fn bv_to_cut_point_distance(&self, axis: &A, bounding_volume: &B, cut_point: &P) -> Self::Dist {
        (self)(axis, bounding_volume, cut_point)
    }
}

pub trait DistanceBVBV<BA, BB> {
    type Dist;

    fn bv_to_bv_distance(&self, bounding_volume_a: &BA, bounding_volume_b: &BB) -> Self::Dist;
}

impl<BA, BB, F, D> DistanceBVBV<BA, BB> for F where F: Fn(&BA, &BB) -> D {
    type Dist = D;

    fn bv_to_bv_distance(&self, bounding_volume_a: &BA, bounding_volume_b: &BB) -> Self::Dist {
        (self)(bounding_volume_a, bounding_volume_b)
    }
}

pub struct KdvTree<A, P, B, S> {
    axis: Vec<A>,
    shapes: Vec<S>,
    root: KdvNode<P, B>,
}

impl<A, P, B, S> KdvTree<A, P, B, S>
    where B: BoundingVolume<P>,
{
    pub fn build<IA, II, CMF, BVF, CPF, CBF>(
        axis_it: IA,
        shapes_it: II,
        cmp_p: CMF,
        get_bv: BVF,
        mut get_cp: CPF,
        mut cut_bv: CBF,
    )
        -> Result<KdvTree<A, P, B, S>, CBF::Error>
        where IA: IntoIterator<Item = A>,
              II: IntoIterator<Item = S>,
              CMF: CmpPoints<A, P>,
              BVF: GetBoundingVolume<B, S>,
              CPF: GetCutPoint<A, P>,
              CBF: BoundingVolumesCutter<A, P, B, S>,
    {
        let axis: Vec<_> = axis_it.into_iter().collect();
        let shapes: Vec<_> = shapes_it.into_iter().collect();
        let root_shapes: Vec<_> = shapes
            .iter()
            .enumerate()
            .map(|(i, s)| ShapeFragment {
                bounding_volume: get_bv.bounding_volume(s),
                shape_id: i,
            })
            .collect();
        Ok(KdvTree {
            root: KdvNode::build(0, &axis, &shapes, root_shapes, &cmp_p, &get_bv, &mut get_cp, &mut cut_bv)?,
            axis, shapes,
        })
    }

    pub fn intersects<'t, 's, SN, BN, CMF, BVF, CPF, CBF>(
        &'t self,
        shape: &'s SN,
        cmp_p: CMF,
        get_bv: BVF,
        get_cp: CPF,
        cut_bv: CBF,
    )
        -> IntersectIter<'t, 's, A, P, S, B, SN, BN, CMF, CPF, CBF>
        where CMF: CmpPoints<A, P>,
              BVF: GetBoundingVolume<BN, SN>,
              CPF: GetCutPoint<A, P>,
              CBF: BoundingVolumesCutter<A, P, BN, SN>,
    {
        IntersectIter {
            needle: shape,
            axis: &self.axis,
            shapes: &self.shapes,
            queue: vec![TraverseTask::Explore {
                node: &self.root,
                needle_fragment: get_bv.bounding_volume(shape),
            }],
            cmp_p,
            get_cp,
            cut_bv,
        }
    }

    pub fn nearest<'t, SN, BN, D, CMF, BVF, DPF, DVF>(
        &'t self,
        shape: &SN,
        cmp_p: CMF,
        get_bv: BVF,
        dist_cp: DPF,
        dist_bv: DVF,
    )
        -> NearestIter<'t, A, P, S, B, BN, D, CMF, DPF, DVF>
        where CMF: CmpPoints<A, P>,
              BVF: GetBoundingVolume<BN, SN>,
              DPF: DistanceBVCP<A, P, BN, Dist = D>,
              DVF: DistanceBVBV<B, BN, Dist = D>,
              D: PartialEq + PartialOrd,
    {
        NearestIter {
            needle_bv: get_bv.bounding_volume(shape),
            axis: &self.axis,
            shapes: &self.shapes,
            nodes_queue: BinaryHeap::from_iter(iter::once(NearestNode {
                dist: None,
                node: &self.root,
            })),
            neighbours: BinaryHeap::new(),
            visited: HashSet::new(),
            cmp_p,
            dist_cp,
            dist_bv,
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
    fn build<A, S, CMF, BVF, CPF, CBF>(
        depth: usize,
        axis: &[A],
        shapes: &[S],
        mut node_shapes: Vec<ShapeFragment<B>>,
        cmp_p: &CMF,
        get_bv: &BVF,
        get_cp: &mut CPF,
        cut_bv: &mut CBF,
    )
        -> Result<KdvNode<P, B>, CBF::Error>
        where B: BoundingVolume<P>,
              CMF: CmpPoints<A, P>,
              BVF: GetBoundingVolume<B, S>,
              CPF: GetCutPoint<A, P>,
              CBF: BoundingVolumesCutter<A, P, B, S>,
    {
        // locate cut point for coords
        let cut_axis_index = depth % axis.len();
        let cut_axis = &axis[cut_axis_index];
        let maybe_cut_point = get_cp.cut_point(
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
                let owner = shape_owner(&shapes[shape_id], bounding_volume, cut_axis, &cut_point, cmp_p, cut_bv)?;
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
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes, cmp_p, get_bv, get_cp, cut_bv)?),
                }
            } else if right_shapes.is_empty() {
                KdvNodeChildren::OnlyLeft {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes, cmp_p, get_bv, get_cp, cut_bv)?),
                }
            } else {
                KdvNodeChildren::Both {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    left: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes, cmp_p, get_bv, get_cp, cut_bv)?),
                    right: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes, cmp_p, get_bv, get_cp, cut_bv)?),
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

fn shape_owner<A, P, B, S, CMF, CBF>(shape: &S, fragment: B, cut_axis: &A, cut_point: &P, cmp_p: &CMF, cut_bv: &mut CBF) ->
    Result<ShapeOwner<B>, CBF::Error>
    where B: BoundingVolume<P>,
          CMF: CmpPoints<A, P>,
          CBF: BoundingVolumesCutter<A, P, B, S>,
{
    let min_corner = fragment.min_corner();
    let max_corner = fragment.max_corner();
    Ok(match (cmp_p.cmp_points(&cut_axis, &min_corner, cut_point), cmp_p.cmp_points(&cut_axis, &max_corner, cut_point)) {
        (Ordering::Less, Ordering::Less) | (Ordering::Less, Ordering::Equal) =>
            ShapeOwner::Left(fragment),
        (Ordering::Greater, Ordering::Greater) | (Ordering::Equal, Ordering::Greater) =>
            ShapeOwner::Right(fragment),
        _ => if let Some((left_bvol, right_bvol)) = cut_bv.cut(shape, &fragment, cut_axis, cut_point)? {
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

pub struct IntersectIter<'t, 's, A: 't, P: 't, SS: 't, BS: 't, SN: 's, BN, CMF, CPF, CBF> {
    needle: &'s SN,
    axis: &'t [A],
    shapes: &'t [SS],
    queue: Vec<TraverseTask<'t, P, BS, BN>>,
    cmp_p: CMF,
    get_cp: CPF,
    cut_bv: CBF,
}

#[derive(Clone, PartialEq, Debug)]
pub struct Intersection<'t, SS: 't, BS: 't, BN> {
    pub shape: &'t SS,
    pub shape_fragment: &'t BS,
    pub needle_fragment: BN,
}

impl<'t, 's, A, P, SS, BS, SN, BN, CMF, CPF, CBF> Iterator for IntersectIter<'t, 's, A, P, SS, BS, SN, BN, CMF, CPF, CBF>
    where BS: BoundingVolume<P>,
          BN: BoundingVolume<P> + Clone,
          CMF: CmpPoints<A, P>,
          CPF: GetCutPoint<A, P>,
          CBF: BoundingVolumesCutter<A, P, BN, SN>,
{
    type Item = Result<Intersection<'t, SS, BS, BN>, CBF::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        'outer: while let Some(task) = self.queue.pop() {
            match task {
                TraverseTask::Explore { node, needle_fragment, } => {
                    match node.children {
                        KdvNodeChildren::Missing =>
                            (),
                        KdvNodeChildren::OnlyLeft { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, &self.cmp_p, &mut self.cut_bv) {
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
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, &self.cmp_p, &mut self.cut_bv) {
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
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, &self.cmp_p, &mut self.cut_bv) {
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
                        (self.cmp_p.cmp_points(axis, &needle_min, &shape_max) == Ordering::Greater ||
                         self.cmp_p.cmp_points(axis, &needle_max, &shape_min) == Ordering::Less)
                    });
                    if no_intersection {
                        continue;
                    }
                    let axis_total = self.axis.len();
                    for _ in 0 .. axis_total {
                        let cut_axis = &self.axis[axis_counter % axis_total];
                        axis_counter += 1;
                        let maybe_cut_point = self.get_cp.cut_point(
                            cut_axis,
                            iter::once(needle_fragment.min_corner())
                                .chain(iter::once(needle_fragment.max_corner())),
                        );
                        if let Some(cut_point) = maybe_cut_point {
                            match self.cut_bv.cut(self.needle, &needle_fragment, cut_axis, &cut_point) {
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

struct NearestNode<'t, P: 't, B: 't, D> where D: PartialEq + PartialOrd {
    dist: Option<D>,
    node: &'t KdvNode<P, B>
}

pub struct NearestShape<'t, S: 't, D> {
    pub dist: D,
    pub shape: &'t S,
}

pub struct NearestIter<'t, A: 't, P: 't, S: 't, B: 't, BN, D, CMF, DPF, DVF> where D: PartialEq + PartialOrd {
    needle_bv: BN,
    axis: &'t [A],
    shapes: &'t [S],
    nodes_queue: BinaryHeap<NearestNode<'t, P, B, D>>,
    neighbours: BinaryHeap<NearestShape<'t, S, D>>,
    visited: HashSet<*const S>,
    cmp_p: CMF,
    dist_cp: DPF,
    dist_bv: DVF,
}

impl<'t, A, P, S, B, BN, D, CMF, DPF, DVF> Iterator for NearestIter<'t, A, P, S, B, BN, D, CMF, DPF, DVF>
    where B: BoundingVolume<P>,
          BN: BoundingVolume<P>,
          D: PartialEq + PartialOrd,
          CMF: CmpPoints<A, P>,
          DPF: DistanceBVCP<A, P, BN, Dist = D>,
          DVF: DistanceBVBV<B, BN, Dist = D>,
{
    type Item = NearestShape<'t, S, D>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let direction =
                match (self.neighbours.peek_mut(), self.nodes_queue.peek_mut()) {
                    (None, None) =>
                        return None,
                    (Some(top_shape), None) =>
                        Ok(PeekMut::pop(top_shape)),
                    (None, Some(top_node)) =>
                        Err(PeekMut::pop(top_node)),
                    (Some(top_shape), Some(top_node)) => {
                        let shape_closer = if let Some(ref top_node_dist) = top_node.dist {
                            &top_shape.dist < top_node_dist
                        } else {
                            false
                        };
                        if shape_closer {
                            Ok(PeekMut::pop(top_shape))
                        } else {
                            Err(PeekMut::pop(top_node))
                        }
                    },
                };
            match direction {
                Ok(nearest_shape) => {
                    let key = nearest_shape.shape as *const _;
                    if self.visited.contains(&key) {
                        continue
                    } else {
                        self.visited.insert(key);
                        return Some(nearest_shape);
                    }
                },
                Err(nearest_node) => {
                    let shapes = self.shapes;
                    let dist_cp = &self.dist_cp;
                    let dist_bv = &self.dist_bv;
                    {
                        let needle_bv = &self.needle_bv;
                        self.neighbours.extend(
                            nearest_node.node.shapes.iter()
                                .map(|fragment| NearestShape {
                                    dist: dist_bv.bv_to_bv_distance(&fragment.bounding_volume, needle_bv),
                                    shape: &shapes[fragment.shape_id],
                                })
                        );
                    }

                    let mut cut_bv = |_shape: &(), _fragment: &_, _cut_axis: &_, _cut_point: &_| Ok(None);

                    struct NeedleBv<'a, B: 'a>(&'a B);
                    impl<'a, B, P> BoundingVolume<P> for NeedleBv<'a, B> where B: BoundingVolume<P> {
                        fn min_corner(&self) -> P { self.0.min_corner() }
                        fn max_corner(&self) -> P { self.0.max_corner() }
                    }

                    match nearest_node.node.children {
                        KdvNodeChildren::Missing =>
                            (),
                        KdvNodeChildren::OnlyLeft { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(&(), NeedleBv(&self.needle_bv), cut_axis, &cut_point.point, &self.cmp_p, &mut cut_bv) {
                                Ok(ShapeOwner::Me(NeedleBv(..))) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, }),
                                Ok(ShapeOwner::Left(NeedleBv(..))) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, }),
                                Ok(ShapeOwner::Right(NeedleBv(..))) =>
                                    (),
                                Ok(ShapeOwner::Both { .. }) =>
                                    unreachable!(),
                                Err(()) =>
                                    unreachable!(),
                            }
                        },
                        KdvNodeChildren::OnlyRight { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(&(), NeedleBv(&self.needle_bv), cut_axis, &cut_point.point, &self.cmp_p, &mut cut_bv) {
                                Ok(ShapeOwner::Me(NeedleBv(..))) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, }),
                                Ok(ShapeOwner::Left(NeedleBv(..))) =>
                                    (),
                                Ok(ShapeOwner::Right(NeedleBv(..))) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, }),
                                Ok(ShapeOwner::Both { .. }) =>
                                    unreachable!(),
                                Err(()) =>
                                    unreachable!(),
                            }
                        },
                        KdvNodeChildren::Both { ref cut_point, ref left, ref right, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(&(), NeedleBv(&self.needle_bv), cut_axis, &cut_point.point, &self.cmp_p, &mut cut_bv) {
                                Ok(ShapeOwner::Me(NeedleBv(..))) => {
                                    self.nodes_queue.push(NearestNode { dist: None, node: left, });
                                    self.nodes_queue.push(NearestNode { dist: None, node: right, });
                                },
                                Ok(ShapeOwner::Left(NeedleBv(..))) => {
                                    self.nodes_queue.push(NearestNode { dist: None, node: left, });
                                    self.nodes_queue.push(NearestNode {
                                        dist: Some(dist_cp.bv_to_cut_point_distance(cut_axis, &self.needle_bv, &cut_point.point)),
                                        node: right,
                                    });
                                },
                                Ok(ShapeOwner::Right(NeedleBv(..))) => {
                                    self.nodes_queue.push(NearestNode {
                                        dist: Some(dist_cp.bv_to_cut_point_distance(cut_axis, &self.needle_bv, &cut_point.point)),
                                        node: left,
                                    });
                                    self.nodes_queue.push(NearestNode { dist: None, node: right, });
                                },
                                Ok(ShapeOwner::Both { .. }) =>
                                    unreachable!(),
                                Err(()) =>
                                    unreachable!(),
                            }
                        },
                    }
                },
            }
        }
    }
}

impl<'t, P, B, D> Ord for NearestNode<'t, P, B, D> where D: PartialOrd {
    fn cmp(&self, other: &NearestNode<'t, P, B, D>) -> Ordering {
        match (&self.dist, &other.dist) {
            (&None, &None) =>
                Ordering::Equal,
            (&None, &Some(..)) =>
                Ordering::Greater,
            (&Some(..), &None) =>
                Ordering::Less,
            (&Some(ref da), &Some(ref db)) =>
                if da < db {
                    Ordering::Greater
                } else if da > db {
                    Ordering::Less
                } else {
                    Ordering::Equal
                },
        }
    }
}

impl<'t, P, B, D> PartialOrd for NearestNode<'t, P, B, D> where D: PartialOrd {
    fn partial_cmp(&self, other: &NearestNode<'t, P, B, D>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'t, P, B, D> Eq for NearestNode<'t, P, B, D> where D: PartialEq + PartialOrd { }

impl<'t, P, B, D> PartialEq for NearestNode<'t, P, B, D> where D: PartialEq + PartialOrd {
    fn eq(&self, other: &NearestNode<'t, P, B, D>) -> bool {
        self.dist == other.dist
    }
}

impl<'t, S, D> Ord for NearestShape<'t, S, D> where D: PartialOrd {
    fn cmp(&self, other: &NearestShape<'t, S, D>) -> Ordering {
        if self.dist < other.dist {
            Ordering::Greater
        } else if self.dist > other.dist {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }
}

impl<'t, S, D> PartialOrd for NearestShape<'t, S, D> where D: PartialOrd {
    fn partial_cmp(&self, other: &NearestShape<'t, S, D>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'t, S, D> Eq for NearestShape<'t, S, D> where D: PartialEq + PartialOrd { }

impl<'t, S, D> PartialEq for NearestShape<'t, S, D> where D: PartialEq + PartialOrd {
    fn eq(&self, other: &NearestShape<'t, S, D>) -> bool {
        self.dist == other.dist
    }
}
