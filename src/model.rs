extern crate nalgebra as na;
use self::na::{ Mat3, det };

const EPSILON : f64 = 0.001;

#[derive(Debug)]
enum ConicShape {
    Hyperbola,
    Parabola,
    Ellipse,
    Circle
}

struct ConicSection {
    matrix : Mat3<f64>
}

impl ConicSection {
    fn classify(&self) -> ConicShape {
        let d = det(&self.matrix);

        if d < EPSILON { return ConicShape::Hyperbola; }
        if d > -EPSILON && d < EPSILON { return ConicShape::Parabola; }

        let a = self.matrix.m11;
        let b = self.matrix.m12 * 2.0;
        let c = self.matrix.m22;

        if (a - c).abs() < EPSILON && b.abs() < EPSILON { return ConicShape::Circle; }

        return ConicShape::Ellipse;
    }
}

#[test]
fn it_works2() {
    let m = ConicSection { matrix: na::one() };
    println!("{:?}", m.classify());
}
