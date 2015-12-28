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

#[derive(Debug)]
enum DegenerateConicShape {
    IntersectingLines,
    ParallelLines,
    Point
}

#[derive(Debug)]
struct ConicSection {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
}

impl ConicSection {
    fn new(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> ConicSection {
        ConicSection {
            a: a,
            b: b,
            c: c,
            d: d,
            e: e,
            f: f
        }
    }

    fn matrix(&self) -> Mat3<f64> {
        Mat3::<f64>::new(self.a,       self.b / 2.0, self.d / 2.0,
                         self.b / 2.0, self.c,       self.e / 2.0,
                         self.d / 2.0, self.e / 2.0, self.f       )
    }

    fn classify(&self) -> Result<ConicShape, DegenerateConicShape> {
        let (c, d) = self.classify_det();
        c
    }

    fn classify_det(&self) -> (Result<ConicShape, DegenerateConicShape>, f64) {
        //Calculate determinants of the full 3x3 matrix and the 2x2 top left sub matrix
        let matrix_a33_det = det(&Mat2::<f64>::new(self.a, self.b / 2.0, self.b / 2.0, self.c));
        let det_full = det(&self.matrix());

        //determinate if the conic section is degenerate or not, then determine the precise type of shape
        let class = if det_full > -EPSILON && det_full < EPSILON {
            if matrix_a33_det < EPSILON { Err(DegenerateConicShape::IntersectingLines) }
            else if matrix_a33_det > EPSILON { Err(DegenerateConicShape::Point) }
            else { Err(DegenerateConicShape::ParallelLines) }
        } else {
            if matrix_a33_det < EPSILON { Ok(ConicShape::Hyperbola) }
            else if matrix_a33_det > -EPSILON && matrix_a33_det < EPSILON { Ok(ConicShape::Parabola) }
            else if (self.a - self.c).abs() < EPSILON && self.b.abs() < EPSILON { Ok(ConicShape::Circle) }
            else { Ok(ConicShape::Ellipse) }
        };

        return (class, det_full);
    }

    fn center(&self) -> Vec2<f64> {
        let ac4_bsqr = 4.0 * self.a * self.c - self.b * self.b;

        Vec2::<f64>::new(
            (self.b * self.e - 2.0 * self.c * self.d) / ac4_bsqr,
            (self.d * self.b - 2.0 * self.a * self.e) / ac4_bsqr,
        )
    }

    fn eccentricity(&self) -> Option<f64> {
        match self.classify_det() {
            (Ok(cls), det) => Some(self.eccentricity_fast(cls, det)),
            (Err(e), _) => None
        }

    }

    fn eccentricity_fast(&self, shape: ConicShape, det: f64) -> f64 {
        match shape {
            ConicShape::Parabola => 1.0,
            ConicShape::Circle => 0.0,
            _ => {
                let n = if (det > 0.0) { -1.0 } else { 1.0 };
                let top = 2.0 * ((self.a - self.c).powf(2.0) + self.b.powf(2.0)).sqrt();
                let bot = n * (self.a + self.c) + ( (self.a - self.c).powf(2.0) + self.b.powf(2.0) ).sqrt();
                (top / bot).sqrt()
            }
        }
    }
}

#[test]
fn assert_ellipse_classify_is_ok_ellipse() {
    let conic = ConicSection::new(2.0, -3.0, 4.0, 6.0, -3.0, -4.0);

    match conic.classify() {
        Ok(ConicShape::Ellipse) => assert!(true),
        Ok(a) => { println!("Ok({0:?})", a); assert!(false) },
        Err(a) => { println!("Err({0:?})", a); assert!(false) },
    }
}

#[test]
fn assert_ellipse_center_is_correct_value() {
    let conic = ConicSection::new(2.0, -3.0, 4.0, 6.0, -3.0, -4.0);

    let center = conic.center();
    println!("({0}, {1})", center.x, center.y);

    assert!((center.x - -1.6956).abs() < EPSILON);
    assert!((center.y - -0.2608).abs() < EPSILON);
}

#[test]
fn assert_ellipse_eccentricity_is_correct_value() {
    let conic = ConicSection::new(2.0, -3.0, 4.0, 6.0, -3.0, -4.0);

    let eccentricity = match conic.eccentricity() {
        Some(e) => e,
        None => panic!()
    };
    println!("e {:?}", eccentricity);

    assert!((eccentricity - 0.866442).abs() < EPSILON);
}
