# RYS 核引积分问题深度分析

## 核心矛盾

**现象**：

- 重叠积分(S)和动能积分(T)：完全正确
- 核引积分(V)：有系统误差(0.1-0.2)

**这个现象表明**：

1. 原子坐标获取正确（否则所有积分都会错）
2. OS_mat 使用 Obara-Saika 算法正确
3. Ry_mat 使用 RYS 算法有误

## 关键代码对比

### 重叠积分（OS_mat - 正确）

```rust
pub fn sval() -> f64 {
    let mut sval = 1.0;
    for d in 0..3 {
        let i = self.Ga.lmn[d] as isize;
        let j = self.Gb.lmn[d] as isize;
        sval *= self.sval_d(i, j, d);
    }
    sval
}
```

使用**三维连乘**：分别计算 x、y、z 方向，再相乘。

### 核引积分（Ry_mat - 有误）

```rust
pub fn vval_c(&self, rc: &[f64; 3]) -> f64 {
    ...
    let [xa, ya, za] = self.Ga.xyz;
    let [xb, yb, zb] = self.Gb.xyz;
    let xp = (a * xa + b * xb) / (a + b);
    let yp = (a * ya + b * yb) / (a + b);
    let zp = (a * za + b * zb) / (a + b);
    let rp = [xp, yp, zp];
    let Rpc = maths::calc_dist(&rp, rc);
    let alpha = self.p * Rpc.powi(2);
    ...
}
```

**计算重叠中心到原子核的距离 Rpc**，然后计算 alpha = p \* Rpc^2。

## 关键疑点

### 1. 重叠中心计算

```rust
// Ry_mat使用
let xp = (a * xa + b * xb) / (a + b);

// OS_mat使用
let Xpa = -self.Gb.alp / self.p * Xab;
let Xpb = self.Ga.alp / self.p * Xab;
```

重叠中心的定义应该是：`xp = (a*xa + b*xb) / (a + b)`

这个计算是正确的。

### 2. alpha 参数

```rust
let alpha = self.p * Rpc.powi(2);
```

这里 Rpc 是重叠中心到原子核的距离。如果原子核位置错误，Rpc 就会错误。

但是如果原子坐标正确，Rpc 也应该是正确的。

### 3. 最可能的错误：RYS 根数插值精度

让我检查一下 RYS 根数是如何计算的：

在 roots.rs 中，根数是通过多项式插值获得的：

```rust
pub fn rys_roots(gamma: usize, alpha: f64) -> Vec<[f64; 2]> {
    let x = RYSX[gamma].clone();
    let y = RYSY[gamma].clone();
    maths::lag_intpol(&x, &y, alpha)
}
```

使用 5 点拉格朗日插值。

**可能的问题**：

1. 插值区间是否足够？
2. 对于不同的 alpha 值，插值精度是否足够？

### 4. Hermite 多项式递归

另一个可能的问题是 eval_d 函数的实现。让我对比一下：

OS_mat 的 sval_d：

```rust
if i == 0 && j == 0 {
    return f64::sqrt(PI / self.p) * f64::exp(-self.u * Xab.powi(2));
}
```

Ry_mat 的 eval_d：

```rust
if i == 0 && j == 0 {
    return if t == 0 { Kab } else { 0.0 };
}
```

这里计算 s=0 时的基础值。

**问题可能是 eval_d 函数中递归计算错误**。

## 最有希望的排查方向

### 1. 验证 RYS 根数表

检查不同 gamma 值对应的根数表是否正确。特别是：

- RYSX[1], RYSY[1]：对应 s-s 轨道的根数
- RYSX[2], RYSY[2]：对应 s-p 轨道的根数

### 2. 验证 alpha 值

对于 H2 分子（两个 1s 轨道），计算 alpha 的值：

- 两个 1s 轨道：l=m=n=0，所以角动量=0
- gamma = (0+0)/2 + 1 = 1
- 所以使用 gamma=1 的根数表

alpha = p \* Rpc^2

如果 Rpc 接近 0，alpha 也会很小。

### 3. 验证 Hermite 多项式

检查 hval_d 函数的递归关系是否正确。

## 立即验证建议

1. **打印关键中间值**：

```rust
let xp = (a * xa + b * xb) / (a + b);
println!("重叠中心: ({}, {}, {})", xp, yp, zp);
let Rpc = maths::calc_dist(&rp, rc);
println!("距离原子核: {}", Rpc);
println!("alpha: {}", alpha);
```

2. **验证 RYS 根数**：

```rust
let roost = rys_roots(gamma, alpha);
println!("RYS根数: {:?}", roost);
```

3. **验证每个维度的贡献**：

```rust
for d in 0..3 {
    let i = self.Ga.lmn[d] as isize;
    let j = self.Gb.lmn[d] as isize;
    let xc = rc[d];
    let ival = self.ival_d(i, j, d, xc, s2);
    println!("维度{}: ival = {}", d, ival);
}
```

## 结论

**最可能的错误原因**：

1. RYS 根数插值表可能不完整或不准确
2. Hermite 多项式递归实现可能有错误
3. 最终系数 `2π/p` 可能不正确

建议先检查 RYS 根数表和插值函数的准确性。
