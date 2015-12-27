extern crate nalgebra as na;
use self::na::{ Mat2, Mat3, Vec2, det, eigen_qr };

const EPSILON : f64 = 0.001;

#[derive(Debug)]
enum ConicShape {
    Hyperbola,
    Parabola,
    Ellipse,
    Circle
}

enum DegenerateConicShape {
    IntersectingLines,
    ParallelLines,
    Point
}

struct ConicSection {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
}

impl ConicSection {
    fn matrix(&self) -> Mat3<f64> {
        Mat3::<f64>::new(self.a,       self.b / 2.0, self.d / 2.0,
                         self.b / 2.0, self.c,       self.e / 2.0,
                         self.d / 2.0, self.e / 2.0, self.f)
    }

    fn classify(&self) -> Result<ConicShape, DegenerateConicShape> {
        //Calculate determinants of the full 3x3 matrix and the 2x2 top left sub matrix
        let matrix_a33_det = det(&Mat2::<f64>::new(self.a, self.b / 2.0, self.b / 2.0, self.c));
        let det_full = det(&self.matrix());

        //determinate if the conic section is degenerate or not
        if det_full > -EPSILON && det_full < EPSILON {
            if matrix_a33_det < EPSILON { return Err(DegenerateConicShape::IntersectingLines); }
            if matrix_a33_det > -EPSILON && matrix_a33_det < EPSILON { return Err(DegenerateConicShape::ParallelLines); }
            return Err(DegenerateConicShape::Point);
        }

        //Not degenerate, measure type
        if matrix_a33_det < EPSILON { return Ok(ConicShape::Hyperbola); }
        if matrix_a33_det > -EPSILON && matrix_a33_det < EPSILON { return Ok(ConicShape::Parabola); }
        if (self.a - self.c).abs() < EPSILON && self.b.abs() < EPSILON { return Ok(ConicShape::Circle); }
        return Ok(ConicShape::Ellipse);
    }

    fn center(&self) -> Vec2<f64> {
        let ac4_bsqr = 4.0 * self.a * self.c - self.b * self.b;

        Vec2::<f64>::new(
            (self.b * self.e - 2.0 * self.c * self.d) / ac4_bsqr,
            (self.d * self.b - 2.0 * self.a * self.e) / ac4_bsqr,
        )
    }
}

#[test]
fn it_works2() {
    println!("{:?}", 1.0 / 0.0)
}
