use crate::base::Mole;
use crate::maths;
use crate::moleprop;
use crate::orbtprop;

use ndarray::Array2;

pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }
}
