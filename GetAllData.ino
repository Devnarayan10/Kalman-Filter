#include <Wire.h>
#include <MPU6050_tockn.h>
#undef F
MPU6050 mpu6050(Wire);

// const int MPU_addr = 0x68; // I2C address of the MPU-6050

int16_t AcX, AcY, AcZ, Tmp, GyX, GyY, GyZ;

struct Vector3d {
    double v[3];

    Vector3d() : v{0.0, 0.0, 0.0} {}

    Vector3d(double x, double y, double z) : v{x, y, z} {}

    double& operator[](int i) { return v[i]; }
    const double& operator[](int i) const { return v[i]; }

    Vector3d operator+(const Vector3d& other) const {
        return Vector3d(v[0] + other.v[0], v[1] + other.v[1], v[2] + other.v[2]);
    }

    Vector3d operator-(const Vector3d& other) const {
        return Vector3d(v[0] - other.v[0], v[1] - other.v[1], v[2] - other.v[2]);
    }

    Vector3d operator*(double scalar) const {
        return Vector3d(v[0] * scalar, v[1] * scalar, v[2] * scalar);
    }

    double dot(const Vector3d& other) const {
        return v[0] * other.v[0] + v[1] * other.v[1] + v[2] * other.v[2];
    }
};

struct Matrix3x3 {
    double m[3][3];

    Matrix3x3() {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i][j] = 0.0;
    }

    Matrix3x3(double a00, double a01, double a02,
              double a10, double a11, double a12,
              double a20, double a21, double a22) {
        m[0][0] = a00; m[0][1] = a01; m[0][2] = a02;
        m[1][0] = a10; m[1][1] = a11; m[1][2] = a12;
        m[2][0] = a20; m[2][1] = a21; m[2][2] = a22;
    }

    static Matrix3x3 identity() {
        return Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    }

    Matrix3x3 operator*(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][0] * other.m[0][j] +
                                 m[i][1] * other.m[1][j] +
                                 m[i][2] * other.m[2][j];
            }
        }
        return result;
    }

    Matrix3x3 operator+(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][j] + other.m[i][j];
            }
        }
        return result;
    }

    Matrix3x3 operator-(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][j] - other.m[i][j];
            }
        }
        return result;
    }

    Matrix3x3 operator*(double scalar) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][j] * scalar;
            }
        }
        return result;
    }

        Vector3d operator*(const Vector3d& v) const {
        return Vector3d{
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
        };
    }

    Matrix3x3 transpose() const {
        return Matrix3x3(
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2]
        );
    }

    Matrix3x3 inverse() const {
        Matrix3x3 result = Matrix3x3::identity();
        Matrix3x3 temp = *this;

        for (int i = 0; i < 3; ++i) {
            double maxVal = fabs(temp.m[i][i]);
            int maxRow = i;
            for (int k = i + 1; k < 3; ++k) {
                if (fabs(temp.m[k][i]) > maxVal) {
                    maxVal = fabs(temp.m[k][i]);
                    maxRow = k;
                }
            }

            if (maxRow != i) {
                // Swap rows in both matrices
                for (int j = 0; j < 3; ++j) {
                    double tmp = temp.m[i][j];
                    temp.m[i][j] = temp.m[maxRow][j];
                    temp.m[maxRow][j] = tmp;

                    tmp = result.m[i][j];
                    result.m[i][j] = result.m[maxRow][j];
                    result.m[maxRow][j] = tmp;
                }
            }

            double pivot = temp.m[i][i];
            for (int j = 0; j < 3; ++j) {
                temp.m[i][j] /= pivot;
                result.m[i][j] /= pivot;
            }

            for (int k = 0; k < 3; ++k) {
                if (k != i) {
                    double factor = temp.m[k][i];
                    for (int j = 0; j < 3; ++j) {
                        temp.m[k][j] -= factor * temp.m[i][j];
                        result.m[k][j] -= factor * result.m[i][j];
                    }
                }
            }
        }
        return result;
    }
};

class KalmanFilter {
public:
    KalmanFilter(const Matrix3x3& f, const Matrix3x3& h, const Matrix3x3& q, const Matrix3x3& r)
        : x(Vector3d()), P(Matrix3x3::identity()), F(f), H(h), Q(q), R(r) {}

    void predict() {
        x = F * x;
        P = F * P * F.transpose() + Q;
    }

    void update(const Vector3d& z) {
        Matrix3x3 S = H * P * H.transpose() + R;
        Matrix3x3 K = P * H.transpose() * S.inverse();
        x = x + K * (z - H * x);
        P = (Matrix3x3::identity() - K * H) * P;
    }

    Vector3d getState() const {
        return x;
    }

private:
    Vector3d x;
    Matrix3x3 P, F, H, Q, R;
};
    Matrix3x3 F = Matrix3x3::identity();
    Matrix3x3 H = Matrix3x3::identity();
    Matrix3x3 Q = Matrix3x3::identity()*0.005;  // No scalar multiplication, use identity
    Matrix3x3 R = Matrix3x3::identity()*0.05;  // No scalar multiplication, use identity
    KalmanFilter kf(F, H, Q, R);
void setup() {
    Serial.begin(9600);
    Wire.begin();
    mpu6050.begin();
    mpu6050.calcGyroOffsets(false);
}
long timer = 0;
void loop() {
    mpu6050.update();
    if(millis() - timer > 1000){
      AcX = mpu6050.getAngleX();
      AcY = mpu6050.getAngleY();
      AcZ = mpu6050.getAngleZ();
      


      kf.predict();
      Vector3d z(AcX, AcY, AcZ);
      kf.update(z);

      Vector3d state = kf.getState();

      // Printing accelerometer data and Kalman filter state
      // Serial.print(","); 
      Serial.print(AcX);
      Serial.print(","); 
      Serial.print(AcY);
      Serial.print(","); 
      Serial.print(AcZ);
      Serial.print(","); 
      Serial.print(state[0]); 
      Serial.print(","); 
      Serial.print(state[1]); 
      Serial.print(","); 
      Serial.println(state[2]);
      //timer = millis();
    }
    delay(500);
}
