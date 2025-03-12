// Part XV: One-dimensional Kalman filter (sensor fusion using gyroscope and accelerometer for angle estimation)
/*
The contents of this code and instructions are the intellectual property of Carbon Aeronautics. 
The text and figures in this code and instructions are licensed under a Creative Commons Attribution - Noncommercial - ShareAlike 4.0 International Public Licence. 
This license lets you remix, adapt, and build upon your work non-commercially, as long as you credit Carbon Aeronautics 
(but not in any way that suggests that we endorse you or your use of the work) and license your new creations under the identical terms.
This code and instruction is provided "As Is” without any further warranty. Neither Carbon Aeronautics or the author has any liability to any person or entity 
with respect to any loss or damage caused or declared to be caused directly or indirectly by the instructions contained in this code or by 
the software and hardware described in it. As Carbon Aeronautics has no control over the use, setup, assembly, modification or misuse of the hardware, 
software and information described in this manual, no liability shall be assumed nor accepted for any resulting damage or injury. 
By the act of copying, use, setup or assembly, the user accepts all resulting liability.

1.0  29 December 2022 -  initial release
*/
#include <Arduino.h>
#include <Wire.h>
#include <BasicLinearAlgebra.h>

using namespace BLA;

template <int Dim, typename DType = float>
struct DiagonalMatrix : public MatrixBase<DiagonalMatrix<Dim>, Dim, Dim, DType>
{
    Matrix<Dim, 1, DType> diagonal;

    // For const matrices (ones whose elements can't be modified) you just need to implement this function:
    DType operator()(int row, int col) const
    {
        // If it's on the diagonal and it's not larger than the matrix dimensions then return the element
        if (row == col)
            return diagonal(row);
        else
            // Otherwise return zero
            return 0.0f;
    }

    // If you want to declare a matrix whose elements can be modified then you'll need to define this function:
    // DType& operator()(int row, int col)
};

float d2r = 3.14159265359/180.0, r2d = 180.0/3.14159265359;
float m2mm = 1000.0, mm2m = 0.001;
float RateRoll, RatePitch, RateYaw;
float RateCalibrationRoll, RateCalibrationPitch, RateCalibrationYaw;
int RateCalibrationNumber;
int RateCalibrationTotalNumber = 1000;
float g = 9.80665;
float AccX, AccY, AccZ;
float AccxOffset = 0.0, AccyOffset = 0.0, AcczOffset = 0.0;
float AccXCompensated = 0.0, AccYCompensated = 0.0, AccZCompensated = 0.0;
float PosX = 0.0, PosY = 0.0, PosZ = 0.0;
float VelX = 0.0, VelY = 0.0, VelZ = 0.0;
float AngleRoll, AnglePitch;
float AngleRollIntegration = 0.0, AnglePitchIntegration = 0.0, AngleYawIntegration = 0.0;
uint32_t LoopTimer;
float Ts = 0.004; // [s]
uint32_t Ts_ms = (int) (Ts * 1000); // [ms]
uint32_t Ts_us = (int) (Ts * 1000000); // [us]
float PNstd = 4.0 * pow(10.0, 0.0);
float MNstd = 2.5 * pow(10.0, -1.0);
float KalmanAngleRoll=0, KalmanUncertaintyAngleRoll=2*2;
float KalmanAnglePitch=0, KalmanUncertaintyAnglePitch=2*2;
float Kalman1DOutput[]={0,0};

void imu_signals(void) {
  Wire.beginTransmission(0x68);
  Wire.write(0x1A);
  Wire.write(0x05);
  Wire.endTransmission();
  Wire.beginTransmission(0x68);
  Wire.write(0x1C);
  Wire.write(0x10);
  Wire.endTransmission();
  Wire.beginTransmission(0x68);
  Wire.write(0x3B);
  Wire.endTransmission(); 
  Wire.requestFrom(0x68,6);
  int16_t AccXLSB = Wire.read() << 8 | Wire.read();
  int16_t AccYLSB = Wire.read() << 8 | Wire.read();
  int16_t AccZLSB = Wire.read() << 8 | Wire.read();
  Wire.beginTransmission(0x68);
  Wire.write(0x1B); 
  Wire.write(0x8);
  Wire.endTransmission();     
  Wire.beginTransmission(0x68);
  Wire.write(0x43);
  Wire.endTransmission();
  Wire.requestFrom(0x68,6);

  int16_t GyroX=Wire.read()<<8 | Wire.read();
  int16_t GyroY=Wire.read()<<8 | Wire.read();
  int16_t GyroZ=Wire.read()<<8 | Wire.read();

  // RateRoll = (float) GyroX * (2.0 / 65.5) * d2r;
  // RatePitch = (float) GyroY * (2.0 / 65.5) * d2r;
  // RateYaw = (float) GyroZ * (2.0 / 65.5) * d2r;

  RateRoll = (float) GyroX * (2.0 / 65.5) * d2r;
  RatePitch = (float) GyroY * (2.0 / 65.5) * d2r;
  RateYaw = (float) GyroZ * (2.0 / 65.5) * d2r;

  // AccxOffset = -0.02;
  // AccyOffset = 0.01;
  // AcczOffset = 0.02;
  // AccX = (float) AccXLSB/4096 + AccxOffset;
  // AccY = (float) AccYLSB/4096 + AccyOffset;
  // AccZ = (float) AccZLSB/4096 + AcczOffset;
  // AccX *= g;
  // AccY *= g;
  // AccZ *= g;


  AccxOffset = -0.13665;
  AccyOffset = 0.08335;
  AcczOffset = 0.23665;
  AccX = (float) AccXLSB * g / 4096 + AccxOffset;
  AccY = (float) AccYLSB * g / 4096 + AccyOffset;
  AccZ = (float) AccZLSB * g / 4096 + AcczOffset;

  AngleRoll = atan(AccY/sqrt(AccX*AccX + AccZ*AccZ));
  AnglePitch = -atan(AccX/sqrt(AccY*AccY + AccZ*AccZ));
}

void remove_omega_bias() {
  RateRoll -= RateCalibrationRoll;
  RatePitch -= RateCalibrationPitch;
  RateYaw -= RateCalibrationYaw;
}

void kalman_1d(float KalmanState, float KalmanUncertainty, float KalmanInput, float KalmanMeasurement) {
  KalmanState = KalmanState + Ts*KalmanInput;
  KalmanUncertainty = KalmanUncertainty + Ts*Ts * PNstd*PNstd;
  float KalmanGain = KalmanUncertainty * 1/(1*KalmanUncertainty + MNstd*MNstd);
  KalmanState = KalmanState + KalmanGain * (KalmanMeasurement-KalmanState);
  KalmanUncertainty = (1-KalmanGain) * KalmanUncertainty;
  Kalman1DOutput[0] = KalmanState; 
  Kalman1DOutput[1] = KalmanUncertainty;
}

void compute_thx(float wx) {
  AngleRollIntegration += Ts * wx;
}

void compute_thy(float wy) {
  AnglePitchIntegration += Ts * wy;
}

void compute_thz(float wz) {
  AngleYawIntegration += Ts * wz;
}

void compute_px_vx(float a) {
  VelX += Ts * a;
  PosX += Ts * VelX + 0.5 * Ts * Ts * a;
}

void compute_py_vy(float a) {
  VelY += Ts * a;
  PosY += Ts * VelY + 0.5 * Ts * Ts * a;
}

void compute_pz_vz(float a) {
  VelZ += Ts * a;
  PosZ += Ts * VelZ + 0.5 * Ts * Ts * a;
}

void print_results(float valueX, float valueY, float valueZ) {
  Serial.print(" = [");
  Serial.print(valueX);
  Serial.print(", ");
  Serial.print(valueY);
  Serial.print(", ");
  Serial.print(valueZ);
  Serial.print("]; ");
}


// Rotation about Z axis
BLA::Matrix<3, 3> Rotz(float thz) {
  BLA::Matrix<3, 3> mat;
  // [ cos(thz)  -sin(thz)   0 ]
  // [ sin(thz)   cos(thz)   0 ]
  // [    0          0       1 ]
  mat(0,0) = cos(thz);
  mat(0,1) = -sin(thz);
  mat(0,2) = 0;
  mat(1,0) = sin(thz);
  mat(1,1) = cos(thz);
  mat(1,2) = 0;
  mat(2,0) = 0;
  mat(2,1) = 0;
  mat(2,2) = 1;
  return mat;
}

// Rotation about Y axis
BLA::Matrix<3, 3> Roty(float thy) {
  BLA::Matrix<3, 3> mat;
  // [  cos(thy)   0   sin(thy) ]
  // [     0       1      0     ]
  // [ -sin(thy)   0   cos(thy) ]
  mat(0,0) = cos(thy);
  mat(0,1) = 0;
  mat(0,2) = sin(thy);
  mat(1,0) = 0;
  mat(1,1) = 1;
  mat(1,2) = 0;
  mat(2,0) = -sin(thy);
  mat(2,1) = 0;
  mat(2,2) = cos(thy);
  return mat;
}

// Rotation about X axis
BLA::Matrix<3, 3> Rotx(float thx) {
  BLA::Matrix<3, 3> mat;
  // [ 1      0         0      ]
  // [ 0   cos(thx)  -sin(thx) ]
  // [ 0   sin(thx)   cos(thx) ]
  mat(0,0) = 1;
  mat(0,1) = 0;
  mat(0,2) = 0;
  mat(1,0) = 0;
  mat(1,1) = cos(thx);
  mat(1,2) = -sin(thx);
  mat(2,0) = 0;
  mat(2,1) = sin(thx);
  mat(2,2) = cos(thx);
  return mat;
}

// Rotation from ZYX Euler angles
BLA::Matrix<3, 3> Rotzyx(float thz, float thy, float thx) {
  return Rotz(thz) * Roty(thy) * Rotx(thx);
}


void setup() {
  Serial.begin(57600);
  pinMode(13, OUTPUT);
  digitalWrite(13, HIGH);
  Wire.setClock(400000);
  Wire.begin();
  delay(250);
  Wire.beginTransmission(0x68); 
  Wire.write(0x6B);
  Wire.write(0x00);
  Wire.endTransmission();
  for (RateCalibrationNumber = 0; RateCalibrationNumber < RateCalibrationTotalNumber; RateCalibrationNumber ++) {
    imu_signals();
    RateCalibrationRoll += RateRoll;
    RateCalibrationPitch += RatePitch;
    RateCalibrationYaw += RateYaw;
    Serial.print("[MPU6050 Gyro Calibration] ");
    Serial.print( (float) RateCalibrationNumber / RateCalibrationTotalNumber * 100);
    Serial.println("%");
    delay(1);
  }
  RateCalibrationRoll /= RateCalibrationTotalNumber;
  RateCalibrationPitch /= RateCalibrationTotalNumber;
  RateCalibrationYaw /= RateCalibrationTotalNumber;
  LoopTimer=micros();
}

void loop() {
  // Get IMU signals
  imu_signals();
  remove_omega_bias();

  // Orientation estimation
  kalman_1d(KalmanAngleRoll, KalmanUncertaintyAngleRoll, RateRoll, AngleRoll);
  KalmanAngleRoll=Kalman1DOutput[0]; 
  KalmanUncertaintyAngleRoll=Kalman1DOutput[1];

  kalman_1d(KalmanAnglePitch, KalmanUncertaintyAnglePitch, RatePitch, AnglePitch);
  KalmanAnglePitch=Kalman1DOutput[0]; 
  KalmanUncertaintyAnglePitch=Kalman1DOutput[1];

  compute_thx(RateRoll);
  compute_thy(RatePitch);
  compute_thz(RateYaw);

  // Position estimation
  // Acc requires orientation estimation to do gravity compensation!!! (e.g. use rotation matrix)
  // Calculate the gravity vector in the sensor frame using rotation matrix from ZYX Euler angles
  BLA::Matrix<3, 3> R_w_s = Rotzyx(AngleYawIntegration, KalmanAnglePitch, KalmanAngleRoll);
  BLA::Matrix<3, 1> g_w(0, 0, -g);
  BLA::Matrix<3, 1> g_s = (~R_w_s) * g_w; // transpose of R_w_s
  BLA::Matrix<3, 1> a_s(AccX, AccY, AccZ);
  BLA::Matrix<3, 1> a_s_compensated = a_s + g_s;

  compute_px_vx(a_s_compensated(0));
  compute_py_vy(a_s_compensated(1));
  compute_pz_vz(a_s_compensated(2));

  // Print results
  // Serial.print("[MPU6050 position and orientation estimation] ");

  // Serial.print("a_s [m/s^2]");
  // print_results(a_s(0), a_s(1), a_s(2));

  // Serial.print("g_s [m/s^2]");
  // print_results(g_s(0), g_s(1), g_s(2));

  // Serial.print("a_s_compensated [m/s^2]");
  // print_results(a_s_compensated(0), a_s_compensated(1), a_s_compensated(2));

  // Serial.print("Vel [mm/s] = [");
  // print_results(VelX * m2mm, VelY * m2mm, VelZ * m2mm);

  // Serial.print("Pos [mm] = [");
  // print_results(PosX * m2mm, PosY * m2mm, PosZ * m2mm);

  // Serial.print("Omega [°/s]");
  // print_results(RateRoll, RatePitch, RateYaw);

  Serial.print("Orientation [°]");
  // print_results(KalmanAngleRoll * r2d, KalmanAnglePitch * r2d, AngleYawIntegration * r2d);
  print_results(AngleRollIntegration * r2d, AnglePitchIntegration * r2d, AngleYawIntegration * r2d);

  Serial.println("");

  while (micros() - LoopTimer < Ts_us);
  LoopTimer=micros();
}