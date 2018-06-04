#[cfg(test)] extern crate rand;

use std::cmp::Ordering;
use std::iter::{self, FromIterator};
use std::collections::{BinaryHeap, binary_heap::PeekMut};

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

    fn bv_to_cut_point_distance(&mut self, axis: &A, bounding_volume: &B, cut_point: &P) -> Self::Dist;
}

impl<A, P, B, F, D> DistanceBVCP<A, P, B> for F where F: FnMut(&A, &B, &P) -> D {
    type Dist = D;

    fn bv_to_cut_point_distance(&mut self, axis: &A, bounding_volume: &B, cut_point: &P) -> Self::Dist {
        (self)(axis, bounding_volume, cut_point)
    }
}

pub trait DistanceBVBV<BA, BB> {
    type Dist;

    fn bv_to_bv_distance(&mut self, bounding_volume_a: &BA, bounding_volume_b: &BB) -> Self::Dist;
}

impl<BA, BB, F, D> DistanceBVBV<BA, BB> for F where F: FnMut(&BA, &BB) -> D {
    type Dist = D;

    fn bv_to_bv_distance(&mut self, bounding_volume_a: &BA, bounding_volume_b: &BB) -> Self::Dist {
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

        struct Op<P, B> {
            node_shapes: Vec<ShapeFragment<B>>,
            instruction: Instruction<P>,
        }

        enum Instruction<P> {
            MakeNode { depth: usize, },
            AssembleOnlyLeft { cut_point: CutPoint<P>, },
            AssembleOnlyRight { cut_point: CutPoint<P>, },
            AssembleBoth { cut_point: CutPoint<P>, },
        }

        let mut ops_stack = vec![Op { node_shapes: root_shapes, instruction: Instruction::MakeNode { depth: 0, }, }];
        let mut ret_stack: Vec<KdvNode<P, B>> = vec![];
        while let Some(op) = ops_stack.pop() {
            match op {
                Op { mut node_shapes, instruction: Instruction::MakeNode { depth, }, } => {
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
                            let owner = shape_owner(&shapes[shape_id], bounding_volume, cut_axis, &cut_point, &cmp_p, &mut cut_bv)?;
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
                        if left_shapes.is_empty() && right_shapes.is_empty() {
                            ret_stack.push(KdvNode {
                                shapes: node_shapes,
                                children: KdvNodeChildren::Missing,
                            });
                        } else if left_shapes.is_empty() {
                            ops_stack.push(Op {
                                node_shapes,
                                instruction: Instruction::AssembleOnlyRight {
                                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                                },
                            });
                            ops_stack.push(Op {
                                node_shapes: right_shapes,
                                instruction: Instruction::MakeNode { depth: depth + 1, },
                            });
                        } else if right_shapes.is_empty() {
                            ops_stack.push(Op {
                                node_shapes,
                                instruction: Instruction::AssembleOnlyLeft {
                                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                                },
                            });
                            ops_stack.push(Op {
                                node_shapes: left_shapes,
                                instruction: Instruction::MakeNode { depth: depth + 1, },
                            });
                        } else {
                            ops_stack.push(Op {
                                node_shapes,
                                instruction: Instruction::AssembleBoth {
                                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                                },
                            });
                            ops_stack.push(Op {
                                node_shapes: left_shapes,
                                instruction: Instruction::MakeNode { depth: depth + 1, },
                            });
                            ops_stack.push(Op {
                                node_shapes: right_shapes,
                                instruction: Instruction::MakeNode { depth: depth + 1, },
                            });
                        }
                    } else {
                        // no cut point choosen, keep all shapes in current node
                        ret_stack.push(KdvNode {
                            shapes: node_shapes,
                            children: KdvNodeChildren::Missing,
                        });
                    }
                },
                Op { node_shapes, instruction: Instruction::AssembleOnlyLeft { cut_point, }, } => {
                    let child = ret_stack.pop().map(Box::new).unwrap_or_else(|| unreachable!());
                    ret_stack.push(KdvNode {
                        shapes: node_shapes,
                        children: KdvNodeChildren::OnlyLeft { cut_point, child },
                    });
                },
                Op { node_shapes, instruction: Instruction::AssembleOnlyRight { cut_point, }, } => {
                    let child = ret_stack.pop().map(Box::new).unwrap_or_else(|| unreachable!());
                    ret_stack.push(KdvNode {
                        shapes: node_shapes,
                        children: KdvNodeChildren::OnlyRight { cut_point, child, },
                    });
                },
                Op { node_shapes, instruction: Instruction::AssembleBoth { cut_point, }, } => {
                    let left = ret_stack.pop().map(Box::new).unwrap_or_else(|| unreachable!());
                    let right = ret_stack.pop().map(Box::new).unwrap_or_else(|| unreachable!());
                    ret_stack.push(KdvNode {
                        shapes: node_shapes,
                        children: KdvNodeChildren::Both { cut_point, left, right },
                    });
                },
            }
        }

        let root = ret_stack.pop().unwrap_or_else(|| unreachable!());
        assert!(ret_stack.is_empty());
        Ok(KdvTree { root, axis, shapes, })
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

    pub fn nearest<'t, 's, SN, BN, D, CMF, BVF, CBF, DPF, DVF>(
        &'t self,
        shape: &'s SN,
        cmp_p: CMF,
        get_bv: BVF,
        cut_bv: CBF,
        dist_cp: DPF,
        dist_bv: DVF,
    )
        -> NearestIter<'t, 's, A, P, S, B, SN, BN, D, CMF, CBF, DPF, DVF>
        where CMF: CmpPoints<A, P>,
              BVF: GetBoundingVolume<BN, SN>,
              CBF: BoundingVolumesCutter<A, P, BN, SN>,
              DPF: DistanceBVCP<A, P, BN, Dist = D>,
              DVF: DistanceBVBV<B, BN, Dist = D>,
              D: PartialEq + PartialOrd,
    {
        NearestIter {
            shape,
            axis: &self.axis,
            shapes: &self.shapes,
            nodes_queue: BinaryHeap::from_iter(iter::once(NearestNode {
                dist: None,
                node: &self.root,
                needle_bv: get_bv.bounding_volume(shape),
            })),
            neighbours: BinaryHeap::new(),
            cmp_p,
            cut_bv,
            dist_cp,
            dist_bv,
        }
    }

    pub fn iter<'t>(&'t self) -> Iter<'t, P, B, S> {
        Iter {
            shapes: &self.shapes,
            pending: vec![(0, &self.root)],
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

struct NearestNode<'t, P: 't, B: 't, BN, D> where D: PartialEq + PartialOrd {
    dist: Option<D>,
    node: &'t KdvNode<P, B>,
    needle_bv: BN,
}

#[derive(Clone, Debug)]
pub struct NearestShape<'t, B: 't, S: 't, D> {
    pub dist: D,
    pub shape: &'t S,
    pub shape_fragment: &'t B,
}

pub struct NearestIter<'t, 's, A: 't, P: 't, S: 't, B: 't, SN: 's, BN, D, CMF, CBF, DPF, DVF> where D: PartialEq + PartialOrd {
    shape: &'s SN,
    axis: &'t [A],
    shapes: &'t [S],
    nodes_queue: BinaryHeap<NearestNode<'t, P, B, BN, D>>,
    neighbours: BinaryHeap<NearestShape<'t, B, S, D>>,
    cmp_p: CMF,
    cut_bv: CBF,
    dist_cp: DPF,
    dist_bv: DVF,
}

impl<'t, 's, A, P, S, B, SN, BN, D, CMF, CBF, DPF, DVF> Iterator for NearestIter<'t, 's, A, P, S, B, SN, BN, D, CMF, CBF, DPF, DVF>
    where B: BoundingVolume<P>,
          BN: BoundingVolume<P> + Clone,
          D: PartialEq + PartialOrd,
          CMF: CmpPoints<A, P>,
          CBF: BoundingVolumesCutter<A, P, BN, SN>,
          DPF: DistanceBVCP<A, P, BN, Dist = D>,
          DVF: DistanceBVBV<B, BN, Dist = D>,
{
    type Item = Result<NearestShape<'t, B, S, D>, CBF::Error>;

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
                Ok(nearest_shape) =>
                    return Some(Ok(nearest_shape)),
                Err(nearest_node) => {
                    let shapes = self.shapes;
                    let dist_cp = &mut self.dist_cp;
                    let dist_bv = &mut self.dist_bv;
                    {
                        let needle_bv = &nearest_node.needle_bv;
                        self.neighbours.extend(
                            nearest_node.node.shapes.iter()
                                .map(|fragment| NearestShape {
                                    dist: dist_bv.bv_to_bv_distance(&fragment.bounding_volume, needle_bv),
                                    shape: &shapes[fragment.shape_id],
                                    shape_fragment: &fragment.bounding_volume,
                                })
                        );
                    }
                    match nearest_node.node.children {
                        KdvNodeChildren::Missing =>
                            (),
                        KdvNodeChildren::OnlyLeft { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.shape, nearest_node.needle_bv, cut_axis, &cut_point.point, &self.cmp_p, &mut self.cut_bv) {
                                Ok(ShapeOwner::Me(needle_fragment)) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, needle_bv: needle_fragment, }),
                                Ok(ShapeOwner::Left(needle_fragment)) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, needle_bv: needle_fragment, }),
                                Ok(ShapeOwner::Right(..)) =>
                                    (),
                                Ok(ShapeOwner::Both { left_bvol: needle_fragment, .. }) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, needle_bv: needle_fragment, }),
                                Err(error) => {
                                    self.neighbours.clear();
                                    self.nodes_queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        },
                        KdvNodeChildren::OnlyRight { ref cut_point, ref child, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.shape, nearest_node.needle_bv, cut_axis, &cut_point.point, &self.cmp_p, &mut self.cut_bv) {
                                Ok(ShapeOwner::Me(needle_fragment)) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, needle_bv: needle_fragment, }),
                                Ok(ShapeOwner::Left(..)) =>
                                    (),
                                Ok(ShapeOwner::Right(needle_fragment)) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, needle_bv: needle_fragment, }),
                                Ok(ShapeOwner::Both { right_bvol: needle_fragment, .. }) =>
                                    self.nodes_queue.push(NearestNode { dist: None, node: child, needle_bv: needle_fragment, }),
                                Err(error) => {
                                    self.neighbours.clear();
                                    self.nodes_queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        },
                        KdvNodeChildren::Both { ref cut_point, ref left, ref right, } => {
                            let cut_axis = &self.axis[cut_point.axis_index];
                            match shape_owner(self.shape, nearest_node.needle_bv, cut_axis, &cut_point.point, &self.cmp_p, &mut self.cut_bv) {
                                Ok(ShapeOwner::Me(needle_fragment)) => {
                                    self.nodes_queue.push(NearestNode { dist: None, node: left, needle_bv: needle_fragment.clone(), });
                                    self.nodes_queue.push(NearestNode { dist: None, node: right, needle_bv: needle_fragment, });
                                },
                                Ok(ShapeOwner::Left(needle_fragment)) => {
                                    self.nodes_queue.push(NearestNode { dist: None, node: left, needle_bv: needle_fragment.clone(), });
                                    self.nodes_queue.push(NearestNode {
                                        dist: Some(dist_cp.bv_to_cut_point_distance(cut_axis, &needle_fragment, &cut_point.point)),
                                        node: right,
                                        needle_bv: needle_fragment,
                                    });
                                },
                                Ok(ShapeOwner::Right(needle_fragment)) => {
                                    self.nodes_queue.push(NearestNode {
                                        dist: Some(dist_cp.bv_to_cut_point_distance(cut_axis, &needle_fragment, &cut_point.point)),
                                        node: left,
                                        needle_bv: needle_fragment.clone(),
                                    });
                                    self.nodes_queue.push(NearestNode { dist: None, node: right, needle_bv: needle_fragment, });
                                },
                                Ok(ShapeOwner::Both { left_bvol, right_bvol, }) => {
                                    self.nodes_queue.push(NearestNode { dist: None, node: left, needle_bv: left_bvol, });
                                    self.nodes_queue.push(NearestNode { dist: None, node: right, needle_bv: right_bvol, });
                                },
                                Err(error) => {
                                    self.neighbours.clear();
                                    self.nodes_queue.clear();
                                    return Some(Err(error));
                                },
                            }
                        },
                    }
                },
            }
        }
    }
}

impl<'t, P, B, BN, D> Ord for NearestNode<'t, P, B, BN, D> where D: PartialOrd {
    fn cmp(&self, other: &NearestNode<'t, P, B, BN, D>) -> Ordering {
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

impl<'t, P, B, BN, D> PartialOrd for NearestNode<'t, P, B, BN, D> where D: PartialOrd {
    fn partial_cmp(&self, other: &NearestNode<'t, P, B, BN, D>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'t, P, B, BN, D> Eq for NearestNode<'t, P, B, BN, D> where D: PartialEq + PartialOrd { }

impl<'t, P, B, BN, D> PartialEq for NearestNode<'t, P, B, BN, D> where D: PartialEq + PartialOrd {
    fn eq(&self, other: &NearestNode<'t, P, B, BN, D>) -> bool {
        self.dist == other.dist
    }
}

impl<'t, B, S, D> Ord for NearestShape<'t, B, S, D> where D: PartialOrd {
    fn cmp(&self, other: &NearestShape<'t, B, S, D>) -> Ordering {
        if self.dist < other.dist {
            Ordering::Greater
        } else if self.dist > other.dist {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }
}

impl<'t, B, S, D> PartialOrd for NearestShape<'t, B, S, D> where D: PartialOrd {
    fn partial_cmp(&self, other: &NearestShape<'t, B, S, D>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'t, B, S, D> Eq for NearestShape<'t, B, S, D> where D: PartialEq + PartialOrd { }

impl<'t, B, S, D> PartialEq for NearestShape<'t, B, S, D> where D: PartialEq + PartialOrd {
    fn eq(&self, other: &NearestShape<'t, B, S, D>) -> bool {
        self.dist == other.dist
    }
}

pub struct Iter<'t, P: 't, B: 't, S: 't> {
    shapes: &'t [S],
    pending: Vec<(usize, &'t KdvNode<P, B>)>,
}

impl<'t, P, B, S> Iterator for Iter<'t, P, B, S> {
    type Item = KdvNodeRef<'t, B, S>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((depth, node)) = self.pending.pop() {
            match node.children {
                KdvNodeChildren::Missing =>
                    (),
                KdvNodeChildren::OnlyLeft { ref child, .. } =>
                    self.pending.push((depth + 1, &*child)),
                KdvNodeChildren::OnlyRight { ref child, .. } =>
                    self.pending.push((depth + 1, &*child)),
                KdvNodeChildren::Both { ref left, ref right, .. } => {
                    self.pending.push((depth + 1, &*right));
                    self.pending.push((depth + 1, &*left));
                },
            }
            Some(KdvNodeRef {
                shapes: self.shapes,
                fragments: &node.shapes,
                depth,
            })
        } else {
            None
        }
    }
}

pub struct KdvNodeRef<'t, B: 't, S: 't> {
    shapes: &'t [S],
    fragments: &'t [ShapeFragment<B>],
    depth: usize,
}

impl<'t, B, S> KdvNodeRef<'t, B, S> {
    pub fn depth(&self) -> usize {
        self.depth
    }

    pub fn shapes(self) -> impl Iterator<Item = (&'t S, &'t B)> {
        self.fragments.into_iter()
            .map(move |fragment| (&self.shapes[fragment.shape_id], &fragment.bounding_volume))
    }
}
