创建计算器结构体

```rs
use crate::base::Mole;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }
}
```

性能分析

```sh
samply record -- rswfn.exe
```

子空间中找到一个几何结构的函数，这个函数受到距离越近的原子的影响越大，且和所有原子之间的连线的向量的偏差应该越小
