use std::f64::consts::PI;
use std::fmt;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;

#[derive(Debug, Clone)]
pub struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    pub fn from_arr(vec: [f64; 3]) -> Self {
        Vec3::new(vec[0], vec[1], vec[2])
    }

    pub fn from_vec(vec: &Vec<f64>) -> Self {
        Vec3::new(vec[0], vec[1], vec[2])
    }

    pub fn into_arr(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    pub fn into_vec(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }

    pub fn zero() -> Self {
        Vec3::new(0.0, 0.0, 0.0)
    }

    pub fn rand() -> Self {
        let x: f64 = rand::random();
        let y: f64 = rand::random();
        let z: f64 = rand::random();
        let l = (x * x + y * y + z * z).sqrt();
        Vec3::new(x / l, y / l, z / l)
    }

    pub fn normalize(&self) -> Vec3 {
        let norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        Vec3::new(self.x / norm, self.y / norm, self.z / norm)
    }

    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn cross(&self, vec3: &Vec3) -> Vec3 {
        let [x1, y1, z1] = self.into_arr();
        let [x2, y2, z2] = vec3.into_arr();
        Vec3::new(y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2)
    }

    pub fn dot(&self, vec3: &Vec3) -> f64 {
        self.x * vec3.x + self.y * vec3.y + self.z * vec3.z
    }

    /// 两向量之间的角度 0-1 pi
    pub fn angle(&self, vec3: &Vec3) -> f64 {
        let cos = self.dot(vec3) / (self.length() * vec3.length());
        cos.acos() / PI
    }

    pub fn project_onto(&self, vec3: &Vec3) -> Vec3 {
        vec3 * self.dot(vec3) / vec3.length()
    }

    /// 旋转一个点，self代表point
    pub fn rotate_point(&self, axis: &Vec3, angle: f64) -> Self {
        let axis = axis.normalize();
        let [px, py, pz] = self.into_arr();
        let [dx, dy, dz] = axis.into_arr();
        let cos = angle.cos();
        let sin = angle.sin();
        let x = px * (dx.powi(2) * (1.0 - cos) + cos) + py * (dx * dy * (1.0 - cos) - dz * sin) + pz * (dx * dz * (1.0 - cos) + dy * sin);
        let y = px * (dx * dy * (1.0 - cos) + dz * sin) + py * (dy.powi(2) * (1.0 - cos) + cos) + pz * (dy * dz * (1.0 - cos) - dx * sin);
        let z = px * (dx * dz * (1.0 - cos) - dy * sin) + py * (dy * dz * (1.0 - cos) + dx * sin) + pz * (dz.powi(2) * (1.0 - cos) + cos);
        Vec3::new(x, y, z)
    }

    /// 旋转箭头，self代表箭头的朝向
    pub fn rotate_arrow(&self, cent: &Vec3, axis: &Vec3, angle: f64) -> Self {
        let point = cent + self;
        let point = point.rotate_point(axis, angle);
        point - cent
    }

    pub fn dist_to_plan(&self, plan: &[f64; 4]) -> f64 {
        let [a, b, c, d] = plan;
        (a * self.x + b * self.y + c * self.z + d) / (a * a + b * b + c * c).sqrt()
    }
}

impl Sub<&Vec3> for &Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Sub<Vec3> for Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Sub<&Vec3> for Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Sub<Vec3> for &Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Add<&Vec3> for &Vec3 {
    type Output = Vec3;
    fn add(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Add<Vec3> for Vec3 {
    type Output = Vec3;
    fn add(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Add<&Vec3> for Vec3 {
    type Output = Vec3;
    fn add(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Add<Vec3> for &Vec3 {
    type Output = Vec3;
    fn add(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, scalar: f64) -> Vec3 {
        Vec3::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Mul<f64> for &Vec3 {
    type Output = Vec3;
    fn mul(self, scalar: f64) -> Vec3 {
        Vec3::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, scalar: f64) -> Vec3 {
        Vec3::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
        self.z *= scalar;
    }
}

impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({:>10.4}, {:>10.4}, {:>10.4}):{:>10.4}", self.x, self.y, self.z, self.length())?;
        Ok(())
    }
}
