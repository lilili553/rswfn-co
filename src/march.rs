use crate::base::Vec3;
use crate::datas::march::MARCH;
use crate::maths::calc_dist;
use ndarray::Array4;
use rayon::prelude::*;
use std::collections::HashMap;

pub struct Cube {
    p0: [f64; 3],
    p1: [f64; 3],
}

impl Cube {
    pub fn new(p0: [f64; 3], p1: [f64; 3]) -> Self {
        let [x0, y0, z0] = p0;
        let [x1, y1, z1] = p1;
        let lx = x1 - x0;
        let ly = y1 - y0;
        let lz = z1 - z0;
        assert!(lx >= 0.0 && ly >= 0.0 && lz >= 0.0, "边长不能为负数");
        Self { p0: p0, p1: p1 }
    }

    // 根据要求的数量计算步长，数量可以适当调整，三个方向的步长
    pub fn calc_step(&self, count: usize) -> [f64; 3] {
        let [x0, y0, z0] = self.p0;
        let [x1, y1, z1] = self.p1;
        let lx = x1 - x0;
        let ly = y1 - y0;
        let lz = z1 - z0;
        assert!(lx > 0.0 && ly > 0.0 && lz > 0.0, "边长不能为负数");
        assert_ne!(count, 0, "数量不能为0");
        let step = ((lx * ly * lz) / count as f64).powf(1.0 / 3.0);
        let nx = (lx / step).ceil() as usize;
        let ny = (ly / step).ceil() as usize;
        let nz = (lz / step).ceil() as usize;
        let sx = lx / nx as f64;
        let sy = ly / ny as f64;
        let sz = lz / nz as f64;
        [sx, sy, sz]
    }

    pub fn volum(&self) -> f64 {
        let [x0, y0, z0] = self.p0;
        let [x1, y1, z1] = self.p1;
        let lx = x1 - x0;
        let ly = y1 - y0;
        let lz = z1 - z0;
        lx * ly * lz
    }

    pub fn p0(&self) -> Vec3 {
        Vec3::new(self.p0[0], self.p0[1], self.p0[2])
    }

    pub fn p1(&self) -> Vec3 {
        Vec3::new(self.p1[0], self.p1[1], self.p1[2])
    }

    // 将一个大的cube提取出几个小的cube
    pub fn grids_v1(&self, count: usize) -> ([usize; 3], Vec<[f64; 3]>, f64) {
        let [x0, y0, z0] = self.p0;
        let [x1, y1, z1] = self.p1;
        let lx = x1 - x0;
        let ly = y1 - y0;
        let lz = z1 - z0;
        assert!(lx > 0.0 && ly > 0.0 && lz > 0.0, "边长不能为负数");
        assert_ne!(count, 0, "数量不能为0");
        let step = ((lx * ly * lz) / count as f64).powf(1.0 / 3.0);
        let nx = (lx / step).ceil() as usize;
        let ny = (ly / step).ceil() as usize;
        let nz = (lz / step).ceil() as usize;
        let sx = lx / nx as f64;
        let sy = ly / ny as f64;
        let sz = lz / nz as f64;

        let shape = [nx, ny, nz];
        let mut grids = Vec::with_capacity((nx * ny * nz) as usize);

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    grids.push([x0 + sx * i as f64, y0 + sy * j as f64, z0 + sz * k as f64]);
                }
            }
        }
        let volum = sx * sy * sz;
        (shape, grids, volum)
    }

    pub fn grids_v2(&self, step: f64) -> ([usize; 3], Vec<[f64; 3]>) {
        let [x0, y0, z0] = self.p0;
        let [x1, y1, z1] = self.p1;
        let lx = x1 - x0;
        let ly = y1 - y0;
        let lz = z1 - z0;

        let nx = (lx / step).ceil() as usize + 1;
        let ny = (ly / step).ceil() as usize + 1;
        let nz = (lz / step).ceil() as usize + 1;
        let sx = lx / (nx - 1) as f64;
        let sy = ly / (ny - 1) as f64;
        let sz = lz / (nz - 1) as f64;

        let shape = [nx, ny, nz];
        let mut grids = Vec::with_capacity((nx * ny * nz) as usize);

        for i in 0..nx {
            let x = x0 + sx * i as f64;
            for j in 0..ny {
                let y = y0 + sy * j as f64;
                for k in 0..nz {
                    let z = z0 + sz * k as f64;
                    grids.push([x, y, z]);
                }
            }
        }
        (shape, grids)
    }
}

// 将正的和负的融合到一块
pub fn value_merge(pverts: &Vec<[f64; 3]>, pfaces: &Vec<[usize; 3]>, nverts: &Vec<[f64; 3]>, nfaces: &Vec<[usize; 3]>) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<isize>) {
    let pverts_count = pverts.len();
    let nverts_count = nverts.len();
    let pfaces_count = pfaces.len();
    let nfaces_count = nfaces.len();

    let total_verts = pverts_count + nverts_count;
    let total_faces = pfaces_count + nfaces_count;
    let mut types = vec![0; total_verts]; // 每个顶点的类型，正的为1，负的为-1
    types[0..pverts_count].fill(1);
    types[pverts_count..total_verts].fill(-1);

    let mut verts = Vec::with_capacity(total_verts);
    verts.extend_from_slice(pverts);
    verts.extend_from_slice(nverts);

    let mut faces = Vec::with_capacity(total_faces);
    faces.extend_from_slice(pfaces);

    for i in 0..nfaces_count {
        let new_face = [nfaces[i][0] + pverts_count as usize, nfaces[i][1] + pverts_count as usize, nfaces[i][2] + pverts_count as usize];
        faces.push(new_face);
    }

    (verts, faces, types)
}

pub fn flatten(
    // 将三维的点拉成一维的
    verts: &Vec<[f64; 3]>,
    faces: &Vec<[usize; 3]>,
) -> (
    Vec<f64>,   // 正顶点
    Vec<usize>, // 正面
) {
    let mut verts_flat = vec![0.0; verts.len() * 3];
    let mut faces_flat = vec![0; faces.len() * 3];

    for i in 0..verts.len() {
        verts_flat[i * 3 + 0] = verts[i][0];
        verts_flat[i * 3 + 1] = verts[i][1];
        verts_flat[i * 3 + 2] = verts[i][2];
    }

    for i in 0..faces.len() {
        faces_flat[i * 3 + 0] = faces[i][0];
        faces_flat[i * 3 + 1] = faces[i][1];
        faces_flat[i * 3 + 2] = faces[i][2];
    }

    (verts_flat, faces_flat)
}

pub fn march_cube_v2<F>(
    func: F,
    cube: Cube,
    isov: f64,
    step: f64,
    xyzs: Option<&Vec<[f64; 3]>>,
    limit: Option<f64>, // 默认0.0，所有的都计算
) -> (Vec<[f64; 3]>, Vec<isize>)
where
    F: Fn(&[f64; 3]) -> f64 + Sync + Send,
{
    let limit = match limit {
        Some(limit) => limit,
        None => 0.0,
    };
    let xyzs = match xyzs {
        Some(xyzs) => xyzs,
        None => &Vec::new(),
    };

    let mut verts = Vec::new();
    let mut types = Vec::new();
    // 先将空间分为几个较大的区域，以1埃为间隔，计算小立方体区域是否与等值面有交叉
    // 每个面上以较小的step计算数值，然后根据数值判断是否与等值面有交叉
    let [x0, y0, z0] = cube.p0.map(|e| e as f64);
    let [x1, y1, z1] = cube.p1.map(|e| e as f64);
    let lx = (x1 - x0) as usize;
    let ly = (y1 - y0) as usize;
    let lz = (z1 - z0) as usize;
    log::info!("等值面提取:{}{}{}", lx, ly, lz);
    // 每个面由两个点确定，如何不重复地遍历所有的面，每个面的key可以由面的四个顶点的索引确定

    let div = 1.0_f64.div_euclid(step) as usize; // 商
    let rem = 1.0_f64.rem_euclid(step); // 余数
    let mut rect_xyzs: HashMap<(char, usize, usize, usize), Vec<[f64; 3]>> = HashMap::new();
    // 找出每个面的格点
    // 对x轴的面遍历
    for ix in 0..(lx + 1) {
        for iy in 0..ly {
            for iz in 0..lz {
                // let mut rect = Array2::<f64>::zeros((div, div));
                let mut grids = Vec::new();
                for i in 0..div {
                    for j in 0..div {
                        let x = x0 + ix as f64;
                        let y = y0 + iy as f64 + i as f64 * step + rem / 2.0;
                        let z = z0 + iz as f64 + j as f64 * step + rem / 2.0;
                        grids.push([x, y, z]);
                        // let v = func(&[x, y, z]);
                        // rect[(i, j)] = v;
                    }
                }
                let key = ('x', ix, iy, iz);
                rect_xyzs.insert(key, grids);
            }
        }
    }
    log::info!("x轴遍历完毕");
    // 对y轴的面遍历
    for ix in 0..lx {
        for iy in 0..(ly + 1) {
            for iz in 0..lz {
                // let mut rect = Array2::<f64>::zeros((div, div));
                let mut grids = Vec::new();
                for i in 0..div {
                    for j in 0..div {
                        let x = x0 + ix as f64 + i as f64 * step + rem / 2.0;
                        let y = y0 + iy as f64;
                        let z = z0 + iz as f64 + j as f64 * step + rem / 2.0;
                        // let v = func(&[x, y, z]);
                        // rect[(i, j)] = v;
                        grids.push([x, y, z]);
                    }
                }
                let key = ('y', ix, iy, iz);
                rect_xyzs.insert(key, grids);
            }
        }
    }
    log::info!("y轴遍历完毕");
    // 对z轴的面遍历
    for ix in 0..lx {
        for iy in 0..ly {
            for iz in 0..(lz + 1) {
                // let mut rect = Array2::<f64>::zeros((div, div));
                let mut grids = Vec::new();
                for i in 0..div {
                    for j in 0..div {
                        let x = x0 + ix as f64 + i as f64 * step + rem / 2.0;
                        let y = y0 + iy as f64 + j as f64 * step + rem / 2.0;
                        let z = z0 + iz as f64;
                        // let v = func(&[x, y, z]);
                        // rect[(i, j)] = v;
                        grids.push([x, y, z]);
                    }
                }
                let key = ('z', ix, iy, iz);
                rect_xyzs.insert(key, grids);
            }
        }
    }
    log::info!("z轴遍历完毕");
    // 获取每个矩形的数值

    let rect_vals: HashMap<(char, usize, usize, usize), Vec<f64>> = rect_xyzs
        .into_par_iter()
        .map(|(key, xyzs)| {
            let vals = xyzs.iter().map(|xyz| func(xyz)).collect();
            (key, vals)
        })
        .collect();
    log::info!("矩形数值计算完毕");
    let mut cubes = Vec::with_capacity(lx * ly * lz);

    // 获取所有的小正方体
    for ix in 0..lx {
        for iy in 0..ly {
            for iz in 0..lz {
                let cp0 = [x0 + ix as f64 + 0.0, y0 + iy as f64 + 0.0, z0 + iz as f64 + 0.0];
                let cp1 = [x0 + ix as f64 + 1.0, y0 + iy as f64 + 1.0, z0 + iz as f64 + 1.0];
                let cube = Cube::new(cp0, cp1);
                let key = (ix, iy, iz);
                cubes.push((key, cube));
            }
        }
    }
    log::info!("小正方体生成完毕");
    // 筛选小正方体(判断小正方体和等值面是否相交)
    let joins: Vec<bool> = cubes
        .par_iter()
        .map(|(key, _cube)| {
            let (ix, iy, iz) = key.clone();
            let keys = [('x', ix + 0, iy, iz), ('x', ix + 1, iy, iz), ('y', ix, iy + 0, iz), ('y', ix, iy + 1, iz), ('z', ix, iy, iz + 0), ('z', ix, iy, iz + 1)];
            for key in keys.iter() {
                // let vals = rect_vals.get(key).expect("不存在这样的正方形面");
                let vals = match rect_vals.get(key) {
                    Some(vals) => vals,
                    None => {
                        log::error!("不存在这样的正方形面");
                        // panic!("不存在这样的正方形面");
                        std::process::exit(1);
                    }
                };
                let min = vals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                let max = vals.iter().fold(-f64::INFINITY, |a, &b| a.max(b));
                if min > isov && max < isov {
                    return false;
                }
                if min > -isov && max < -isov {
                    return false;
                }
            }
            return true;
        })
        .collect();
    log::info!("小正方体筛选完毕:{}", cubes.len());

    let par_res: Vec<_> = cubes
        .par_iter()
        .enumerate()
        .filter_map(|(i, (_key, cube))| {
            if joins[i] {
                // 如果是满足条件的正方体
                let (shape, grids) = cube.grids_v2(step);
                // 计算网格点到分子的最近距离
                let dists: Vec<f64> = if limit > 0.0 {
                    grids
                        .iter()
                        .map(|grid| {
                            let mut min_dist = f64::INFINITY;
                            for xyz in xyzs.iter() {
                                let now_dist = calc_dist(xyz, grid);
                                if now_dist < min_dist {
                                    min_dist = now_dist;
                                }
                            }
                            min_dist
                        })
                        .collect()
                } else {
                    vec![0.0; grids.len()]
                };
                let values: Vec<f64> = (0..grids.len())
                    .into_iter()
                    .map(|g| {
                        if limit > 0.0 {
                            // 根据是否有距离限制减小一些计算量
                            if dists[g] > limit { f64::NAN } else { func(&grids[g]) }
                        } else {
                            func(&grids[g])
                        }
                        // func(&grids[g])
                    })
                    .collect();
                let (cverts, ctypes) = march_cube_v1(&shape, &grids, &values, isov);
                Some((cverts, ctypes))
            } else {
                None
            }
        })
        .collect();
    log::info!("小正方体计算完毕");
    for (cverts, ctypes) in par_res {
        verts.extend(cverts);
        types.extend(ctypes);
    }
    (verts, types)
}

pub fn march_cube_v1(
    //等值面提取
    shape: &[usize; 3],
    grids: &Vec<[f64; 3]>,
    values: &Vec<f64>,
    isov: f64,
) -> (Vec<[f64; 3]>, Vec<isize>) {
    let voxel = grids2voxel(shape, grids, values);
    let (pverts, nverts) = voxel2verts(&voxel, isov);
    let mut types = vec![];
    types.extend(vec![1; pverts.len()]);
    types.extend(vec![-1; nverts.len()]);
    let mut verts = vec![];
    verts.extend(pverts);
    verts.extend(nverts);
    (verts, types)
}

// 将格点数据转为体素数据
pub fn grids2voxel(shape: &[usize; 3], grids: &Vec<[f64; 3]>, values: &Vec<f64>) -> Array4<f64> {
    let [nx, ny, nz] = *shape;
    let mut voxel = Array4::zeros((nx as usize, ny as usize, nz as usize, 4));
    let mut n = 0;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let [x, y, z] = grids[n];
                let v = values[n];
                voxel[(i as usize, j as usize, k as usize, 0)] = x;
                voxel[(i as usize, j as usize, k as usize, 1)] = y;
                voxel[(i as usize, j as usize, k as usize, 2)] = z;
                voxel[(i as usize, j as usize, k as usize, 3)] = v;
                n += 1;
            }
        }
    }
    voxel
}

// 将体素数据转为正方形数据
pub fn voxel2cubes(voxel: &Array4<f64>, isov: f64) -> Vec<Cube> {
    let shape = voxel.shape();
    let nx = shape[0];
    let ny = shape[1];
    let nz = shape[2];

    let mut cubes = Vec::new();
    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            for k in 0..nz - 1 {
                // let (xyzs, vals) = get_around(voxel, i, j, k);
                let x0 = voxel[(i, j, k, 0)];
                let y0 = voxel[(i, j, k, 1)];
                let z0 = voxel[(i, j, k, 2)];
                let x1 = voxel[(i + 1, j + 1, k + 1, 0)];
                let y1 = voxel[(i + 1, j + 1, k + 1, 1)];
                let z1 = voxel[(i + 1, j + 1, k + 1, 2)];
                let vals = get_around_v2(voxel, i, j, k);
                // 获取vals的最小值及最小值所在的索引
                let vmin = vals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                let vmax = vals.iter().fold(-f64::INFINITY, |a, &b| a.max(b));

                // 处理isov的等值面
                if vmin < isov && isov < vmax {
                    cubes.push(Cube::new([x0, y0, z0], [x1, y1, z1]));
                }

                // 处理-isov的等值面
                if vmin < -isov && -isov < vmax {
                    cubes.push(Cube::new([x0, y0, z0], [x1, y1, z1]));
                }
            }
        }
    }
    cubes
}

#[rustfmt::skip]
static AROUND2:[[usize;2];19]=[
    [0,1],[1,2],[2,3],[3,0],
    [4,5],[5,6],[6,7],[7,4],
    [0,4],[1,5],[2,6],[4,7],
    [0,2],[0,5],[0,7],
    [6,1],[6,3],[6,4],
    [0,6]
];
fn get_around_v2(voxel: &Array4<f64>, i: usize, j: usize, k: usize) -> Vec<f64> {
    let mut vals = Vec::new();
    for d in 0..8 {
        let [dx, dy, dz] = DIDX[d];
        vals.push(voxel[(i + dx, j + dy, k + dz, 3)]);
    }
    for [d0, d1] in AROUND2 {
        let [dx0, dy0, dz0] = DIDX[d0];
        let [dx1, dy1, dz1] = DIDX[d1];
        let v0 = voxel[(i + dx0, j + dy0, k + dz0, 3)];
        let v1 = voxel[(i + dx1, j + dy1, k + dz1, 3)];
        vals.push((v0 + v1) / 2.0);
    }
    vals
}

// 将体素数据转为格点数据
pub fn voxel2verts(voxel: &Array4<f64>, isov: f64) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let shape = voxel.shape();
    let nx = shape[0];
    let ny = shape[1];
    let nz = shape[2];

    let mut pverts: Vec<[f64; 3]> = Vec::new();
    let mut nverts: Vec<[f64; 3]> = Vec::new();
    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            for k in 0..nz - 1 {
                let (xyzs, vals) = get_around_v1(voxel, i, j, k);
                // let dist = voxel[(i, j, k, 4)]; // 根据当前正方体到分子的做小距离，直接跳过
                let nan = vals.iter().any(|val| val.is_nan());
                if nan {
                    continue;
                }
                // let vals2 = get_around_v2(voxel, i, j, k);
                // 获取vals的最小值及最小值所在的索引
                let vmin = vals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
                let vmax = vals.iter().fold(-f64::INFINITY, |a, &b| a.max(b));

                // 处理isov的等值面
                if vmin < isov && isov < vmax {
                    let key = (0..8).fold(0, |acc, n| acc + if vals[n] > isov { 1 << (7 - n) } else { 0 });
                    pverts.extend(get_boxverts(key, &xyzs, &vals, isov));
                }

                // 处理-isov的等值面
                if vmin < -isov && -isov < vmax {
                    let key_neg = (0..8).fold(0, |acc, n| acc + if vals[n] > -isov { 1 << (7 - n) } else { 0 });
                    let mut verts = get_boxverts(key_neg, &xyzs, &vals, -isov);
                    verts.reverse();
                    nverts.extend(verts);
                }
            }
        }
    }
    (pverts, nverts)
}

static DIDX: [[usize; 3]; 8] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
// 获取一个点周围的点
fn get_around_v1(voxel: &Array4<f64>, i: usize, j: usize, k: usize) -> ([[f64; 3]; 8], [f64; 8]) {
    let mut xyzs = [[0.0; 3]; 8];
    let mut vals = [0.0; 8];

    for n in 0..8 {
        let [di, dj, dk] = DIDX[n];
        let x = voxel[(i + di, j + dj, k + dk, 0)];
        let y = voxel[(i + di, j + dj, k + dk, 1)];
        let z = voxel[(i + di, j + dj, k + dk, 2)];
        let v = voxel[(i + di, j + dj, k + dk, 3)];
        xyzs[n] = [x, y, z];
        vals[n] = v;
    }
    (xyzs, vals)
}

static EDGES: [[isize; 2]; 12] = [[0, 1], [1, 2], [2, 3], [0, 3], [1, 5], [2, 6], [3, 7], [0, 4], [4, 5], [5, 6], [6, 7], [7, 4]];
// 获取等值面与小盒子的交点
fn get_boxverts(key: usize, xyzs: &[[f64; 3]; 8], vals: &[f64; 8], isov: f64) -> Vec<[f64; 3]> {
    let mut verts: Vec<[f64; 3]> = Vec::new();
    let nface = MARCH[key].len() / 3; // 面的数量，每个面3个顶点
    for f in 0..nface {
        let face = MARCH[key][f * 3..(f + 1) * 3].to_vec();
        for p in 0..3 {
            let e = face[p];
            let a = EDGES[e as usize][0];
            let b = EDGES[e as usize][1];
            let pa = xyzs[a as usize];
            let pb = xyzs[b as usize];
            let va = vals[a as usize];
            let vb = vals[b as usize];
            let t = (isov - va) / (vb - va);
            let x = pa[0] + t * (pb[0] - pa[0]);
            let y = pa[1] + t * (pb[1] - pa[1]);
            let z = pa[2] + t * (pb[2] - pa[2]);
            verts.push([x, y, z]);
        }
    }
    verts
}

// 将顶点合并
pub fn merge_verts(old_verts: &Vec<[f64; 3]>) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut verts: Vec<[f64; 3]> = Vec::new();
    let mut faces: Vec<[usize; 3]> = Vec::new();
    let nvert = old_verts.len();
    let mut filts = vec![false; nvert]; // 标记顶点是否已被合并
    let mut fidxs = vec![0; nvert]; // 存储新顶点索引

    for i in 0..nvert {
        if filts[i] {
            continue;
        }
        // 添加当前顶点到新数组，并记录新索引
        let new_idx = verts.len();
        // 计算所有合并顶点的平均位置
        let mut sum = [old_verts[i][0], old_verts[i][1], old_verts[i][2]];
        let mut count = 1;
        filts[i] = true;
        fidxs[i] = new_idx;
        // 遍历所有顶点，合并相近的顶点
        for j in 0..nvert {
            if i == j || filts[j] {
                continue;
            }
            if calc_dist(&old_verts[i], &old_verts[j]) < 0.0001 {
                filts[j] = true;
                fidxs[j] = new_idx;
                sum[0] += old_verts[j][0];
                sum[1] += old_verts[j][1];
                sum[2] += old_verts[j][2];
                count += 1;
            }
        }
        // 使用平均位置作为新顶点位置
        let avg = [sum[0] / count as f64, sum[1] / count as f64, sum[2] / count as f64];
        verts.push(avg);
    }

    // 生成面，确保索引不越界
    for i in (0..nvert).step_by(3) {
        if i + 2 >= nvert {
            break;
        }
        let idx0 = fidxs[i];
        let idx1 = fidxs[i + 1];
        let idx2 = fidxs[i + 2];
        if idx0 != idx1 && idx1 != idx2 && idx0 != idx2 {
            faces.push([idx0 as usize, idx1 as usize, idx2 as usize]);
        }
    }

    (verts, faces)
}
