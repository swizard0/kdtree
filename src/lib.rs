use std::iter;
use std::cmp::Ordering;

#[cfg(test)]
mod tests;

pub trait Axis<P> {
    fn cmp_points(&self, a: &P, b: &P) -> Ordering;
}

pub trait BoundingVolume<P> {
    fn min_corner(&self) -> P;
    fn max_corner(&self) -> P;
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

pub struct KdvTree<A, P, B, S> {
    axis: Vec<A>,
    shapes: Vec<S>,
    root: KdvNode<P, B>,
}

impl<A, P, B, S> KdvTree<A, P, B, S>
    where A: Axis<P>,
          B: BoundingVolume<P>,
{
    pub fn build<IA, II, BVF, CPF, CBF>(
        axis_it: IA,
        shapes_it: II,
        get_bv: BVF,
        mut get_cp: CPF,
        mut cut_bv: CBF,
    )
        -> Result<KdvTree<A, P, B, S>, CBF::Error>
        where IA: IntoIterator<Item = A>,
              II: IntoIterator<Item = S>,
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
            root: KdvNode::build(0, &axis, &shapes, root_shapes, &get_bv, &mut get_cp, &mut cut_bv)?,
            axis, shapes,
        })
    }

    pub fn intersects<'t, 's, SN, BN, BVF, CPF, CBF>(
        &'t self,
        shape: &'s SN,
        get_bv: BVF,
        get_cp: CPF,
        cut_bv: CBF,
    )
        -> IntersectIter<'t, 's, A, P, S, B, SN, BN, CPF, CBF>
        where BVF: GetBoundingVolume<BN, SN>,
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
            get_cp,
            cut_bv,
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
    fn build<A, S, BVF, CPF, CBF>(
        depth: usize,
        axis: &[A],
        shapes: &[S],
        mut node_shapes: Vec<ShapeFragment<B>>,
        get_bv: &BVF,
        get_cp: &mut CPF,
        cut_bv: &mut CBF,
    )
        -> Result<KdvNode<P, B>, CBF::Error>
        where B: BoundingVolume<P>,
              A: Axis<P>,
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
                let owner = shape_owner(&shapes[shape_id], bounding_volume, cut_axis, &cut_point, cut_bv)?;
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
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes, get_bv, get_cp, cut_bv)?),
                }
            } else if right_shapes.is_empty() {
                KdvNodeChildren::OnlyLeft {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    child: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes, get_bv, get_cp, cut_bv)?),
                }
            } else {
                KdvNodeChildren::Both {
                    cut_point: CutPoint { axis_index: cut_axis_index, point: cut_point, },
                    left: Box::new(KdvNode::build(depth + 1, axis, shapes, left_shapes, get_bv, get_cp, cut_bv)?),
                    right: Box::new(KdvNode::build(depth + 1, axis, shapes, right_shapes, get_bv, get_cp, cut_bv)?),
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

fn shape_owner<A, P, B, S, CBF>(shape: &S, fragment: B, cut_axis: &A, cut_point: &P, cut_bv: &mut CBF) ->
    Result<ShapeOwner<B>, CBF::Error>
    where A: Axis<P>,
          B: BoundingVolume<P>,
          CBF: BoundingVolumesCutter<A, P, B, S>,
{
    let min_corner = fragment.min_corner();
    let max_corner = fragment.max_corner();
    Ok(match (cut_axis.cmp_points(&min_corner, cut_point), cut_axis.cmp_points(&max_corner, cut_point)) {
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

pub struct IntersectIter<'t, 's, A: 't, P: 't, SS: 't, BS: 't, SN: 's, BN, CPF, CBF> {
    needle: &'s SN,
    axis: &'t [A],
    shapes: &'t [SS],
    queue: Vec<TraverseTask<'t, P, BS, BN>>,
    get_cp: CPF,
    cut_bv: CBF,
}

#[derive(Clone, PartialEq, Debug)]
pub struct Intersection<'t, SS: 't, BS: 't, BN> {
    pub shape: &'t SS,
    pub shape_fragment: &'t BS,
    pub needle_fragment: BN,
}

impl<'t, 's, A, P, SS, BS, SN, BN, CPF, CBF> Iterator for IntersectIter<'t, 's, A, P, SS, BS, SN, BN, CPF, CBF>
    where A: Axis<P>,
          BS: BoundingVolume<P>,
          BN: BoundingVolume<P> + Clone,
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
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, &mut self.cut_bv) {
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
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, &mut self.cut_bv) {
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
                            match shape_owner(self.needle, needle_fragment.clone(), cut_axis, &cut_point.point, &mut self.cut_bv) {
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
