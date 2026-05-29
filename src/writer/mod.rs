mod cub;
mod dmo;
mod fch;
mod gjf;
mod obj;
mod xyz;

pub use cub::CubWriter;
pub use dmo::DmoWriter;
pub use fch::FchWriter;
pub use gjf::GjfWriter;
pub use obj::ObjWriter;
pub use xyz::XyzWriter;

pub trait Writer {
    fn read_mole(&mut self, mole: &crate::base::Mole);
}
