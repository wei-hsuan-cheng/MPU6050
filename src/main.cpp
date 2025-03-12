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

float RateRoll, RatePitch, RateYaw;
float RateCalibrationRoll, RateCalibrationPitch, RateCalibrationYaw;
int RateCalibrationNumber;
float g = 9.80665;
float AccX, AccY, AccZ;
float AccxOffset, AccyOffset, AcczOffset;
float PosX = 0.0, PosY = 0.0, PosZ = 0.0;
float VelX = 0.0, VelY = 0.0, VelZ = 0.0;
float AngleRoll, AnglePitch;
float AngleYaw = 0.0;
uint32_t LoopTimer;
uint32_t Ts_ms = 4;
float Ts = 0.004;
float PNstd = 4.0 * pow(10, 0.0);
float MNstd = 2.5 * pow(10, -1.0);
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
  RateRoll=(float)GyroX/65.5;
  RatePitch=(float)GyroY/65.5;
  RateYaw=(float)GyroZ/65.5;

  AccxOffset = -0.02;
  AccyOffset = 0.01;
  AcczOffset = 0.02;
  AccX = (float)AccXLSB/4096 + AccxOffset;
  AccY = (float)AccYLSB/4096 + AccyOffset;
  AccZ = (float)AccZLSB/4096 + AcczOffset;
  AccX *= g;
  AccY *= g;
  AccZ *= g;

  AngleRoll = atan(AccY/sqrt(AccX*AccX + AccZ*AccZ)) * 1/(3.142/180);
  AnglePitch = -atan(AccX/sqrt(AccY*AccY + AccZ*AccZ)) * 1/(3.142/180);
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

void compute_thz(float wz) {
  AngleYaw += Ts * wz;
}

void compute_px_vx(float a) {
  VelX += Ts * a;
  PosX += Ts * VelX;
}

void compute_py_vy(float a) {
  VelY += Ts * a;
  PosY += Ts * VelY;
}

void compute_pz_vz(float a) {
  VelZ += Ts * a;
  PosZ += Ts * VelZ;
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

void TestBLA()
{
    // If you've been through the HowToUse example you'll know that you can allocate a Matrix and explicitly specify
    // it's type like so:
    BLA::Matrix<4, 4> mat;

    // And as before it's a good idea to fill the matrix before we use it
    mat.Fill(1);

    // Now let's declare a diagonal matrix. To do that we pass the Diagonal class from above along with whatever
    // template parameters as a template parameter to Matrix, like so:
    DiagonalMatrix<4> diag;

    // If we fill diag we'll get a matrix with all 1's along the diagonal, the identity matrix.
    diag.diagonal.Fill(1);

    // So multiplying it with mat will do nothing:
    Serial.print("still ones: ");
    Serial.println(diag * mat);

    // Diagonal matrices have the handy property of scaling either the rows (premultiplication) or columns
    // (postmultiplication) of a matrix

    // So if we modify the diagonal
    for (int i = 0; i < diag.Rows; i++) diag.diagonal(i) = i + 1;

    // And multiply again, we'll see that the rows have been scaled
    Serial.print("scaled rows: ");
    Serial.print(diag * mat);

    // Point being, if you define a class which serves up something when called upon by the () operator, you can embed
    // it in a matrix and define any kind of behaviour you like. Hopefully that'll let this library support lots more
    // applications while catering to the arduino's limited amount of memory.
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
  for (RateCalibrationNumber=0; RateCalibrationNumber<2000; RateCalibrationNumber ++) {
    imu_signals();
    RateCalibrationRoll+=RateRoll;
    RateCalibrationPitch+=RatePitch;
    RateCalibrationYaw+=RateYaw;
    delay(1);
  }
  RateCalibrationRoll/=2000;
  RateCalibrationPitch/=2000;
  RateCalibrationYaw/=2000;
  LoopTimer=micros();

  // TestBLA();

}

void loop() {
  // Get IMU signals
  imu_signals();

  // Orientation estimation
  RateRoll-=RateCalibrationRoll;
  RatePitch-=RateCalibrationPitch;
  RateYaw-=RateCalibrationYaw;
  
  kalman_1d(KalmanAngleRoll, KalmanUncertaintyAngleRoll, RateRoll, AngleRoll);
  KalmanAngleRoll=Kalman1DOutput[0]; 
  KalmanUncertaintyAngleRoll=Kalman1DOutput[1];

  kalman_1d(KalmanAnglePitch, KalmanUncertaintyAnglePitch, RatePitch, AnglePitch);
  KalmanAnglePitch=Kalman1DOutput[0]; 
  KalmanUncertaintyAnglePitch=Kalman1DOutput[1];

  compute_thz(RateYaw);

  // Position estimation
  // Acc requires orientation estimation to do gravity compensation!!! (e.g. use rotation matrix)
  compute_px_vx(AccX);
  compute_py_vy(AccY);
  compute_pz_vz(AccZ);

  // Print results
  // Serial.print("Acc [m / s^2]");
  // print_results(AccX, AccY, AccZ);

  // Serial.print("Vel [m / s] = [");
  // print_results(VelX, VelY, VelZ);

  // Serial.print("Pos [m] = [");
  // print_results(PosX, PosY, PosZ);

  // Serial.print("Omega [° / s]");
  // print_results(RateRoll, RatePitch, RateYaw);

  Serial.print("Orientation [°]");
  print_results(KalmanAngleRoll, KalmanAnglePitch, AngleYaw);

  Serial.println("");

  while (micros() - LoopTimer < 4000);
  LoopTimer=micros();
}




// // Part XIV: Measure angles with the accelerometer
// /*
// The contents of this code and instructions are the intellectual property of Carbon Aeronautics. 
// The text and figures in this code and instructions are licensed under a Creative Commons Attribution - Noncommercial - ShareAlike 4.0 International Public Licence. 
// This license lets you remix, adapt, and build upon your work non-commercially, as long as you credit Carbon Aeronautics 
// (but not in any way that suggests that we endorse you or your use of the work) and license your new creations under the identical terms.
// This code and instruction is provided "As Is” without any further warranty. Neither Carbon Aeronautics or the author has any liability to any person or entity 
// with respect to any loss or damage caused or declared to be caused directly or indirectly by the instructions contained in this code or by 
// the software and hardware described in it. As Carbon Aeronautics has no control over the use, setup, assembly, modification or misuse of the hardware, 
// software and information described in this manual, no liability shall be assumed nor accepted for any resulting damage or injury. 
// By the act of copying, use, setup or assembly, the user accepts all resulting liability.

// 1.0  29 December 2022 -  initial release
// */
// #include <Arduino.h>
// #include <Wire.h>
// float RateRoll, RatePitch, RateYaw;
// float AccX, AccY, AccZ;
// float AccxOffset, AccyOffset, AcczOffset;
// float AngleRoll, AnglePitch;
// float LoopTimer;
// void gyro_signals(void) {
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1A);
//   Wire.write(0x05);
//   Wire.endTransmission();
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1C);
//   Wire.write(0x10);
//   Wire.endTransmission();
//   Wire.beginTransmission(0x68);
//   Wire.write(0x3B);
//   Wire.endTransmission(); 
//   Wire.requestFrom(0x68,6);
//   int16_t AccXLSB = Wire.read() << 8 | Wire.read();
//   int16_t AccYLSB = Wire.read() << 8 | Wire.read();
//   int16_t AccZLSB = Wire.read() << 8 | Wire.read();
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1B); 
//   Wire.write(0x8);
//   Wire.endTransmission();                                                   
//   Wire.beginTransmission(0x68);
//   Wire.write(0x43);
//   Wire.endTransmission();
//   Wire.requestFrom(0x68,6);
//   int16_t GyroX=Wire.read()<<8 | Wire.read();
//   int16_t GyroY=Wire.read()<<8 | Wire.read();
//   int16_t GyroZ=Wire.read()<<8 | Wire.read();
//   RateRoll=(float)GyroX/65.5;
//   RatePitch=(float)GyroY/65.5;
//   RateYaw=(float)GyroZ/65.5;
//   AccxOffset = -0.02;
//   AccyOffset = 0.01;
//   AcczOffset = 0.02;
//   AccX=(float)AccXLSB/4096 + AccxOffset;
//   AccY=(float)AccYLSB/4096 + AccyOffset;
//   AccZ=(float)AccZLSB/4096 + AcczOffset;
//   AngleRoll=atan(AccY/sqrt(AccX*AccX+AccZ*AccZ))*1/(3.142/180);
//   AnglePitch=-atan(AccX/sqrt(AccY*AccY+AccZ*AccZ))*1/(3.142/180);
// }
// void setup() {
//   Serial.begin(57600);
//   pinMode(13, OUTPUT);
//   digitalWrite(13, HIGH);
//   Wire.setClock(400000);
//   Wire.begin();
//   delay(250);
//   Wire.beginTransmission(0x68); 
//   Wire.write(0x6B);
//   Wire.write(0x00);
//   Wire.endTransmission();
// }
// void loop() {
//   gyro_signals();

//   // Serial.print("[ax, ay, az] [g] = [");
//   // Serial.print(AccX);
//   // Serial.print(", ");
//   // Serial.print(AccY);
//   // Serial.print(", ");
//   // Serial.print(AccZ);
//   // Serial.println("]");

//   Serial.print("[thx, thy] [°] = [");
//   Serial.print(AngleRoll);
//   Serial.print(", ");
//   Serial.print(AnglePitch);
//   Serial.println("]");

//   delay(50);
// }




// // Part V: Gyroscope Calibration
// /*
// The contents of this code and instructions are the intellectual property of Carbon Aeronautics. 
// The text and figures in this code and instructions are licensed under a Creative Commons Attribution - Noncommercial - ShareAlike 4.0 International Public Licence. 
// This license lets you remix, adapt, and build upon your work non-commercially, as long as you credit Carbon Aeronautics 
// (but not in any way that suggests that we endorse you or your use of the work) and license your new creations under the identical terms.
// This code and instruction is provided "As Is” without any further warranty. Neither Carbon Aeronautics or the author has any liability to any person or entity 
// with respect to any loss or damage caused or declared to be caused directly or indirectly by the instructions contained in this code or by 
// the software and hardware described in it. As Carbon Aeronautics has no control over the use, setup, assembly, modification or misuse of the hardware, 
// software and information described in this manual, no liability shall be assumed nor accepted for any resulting damage or injury. 
// By the act of copying, use, setup or assembly, the user accepts all resulting liability.

// 1.0  5 October 2022 -  initial release
// */

// #include <Arduino.h>
// #include <Wire.h>
// float RateRoll, RatePitch, RateYaw;
// float RateCalibrationRoll, RateCalibrationPitch, RateCalibrationYaw;
// int RateCalibrationNumber;
// void gyro_signals(void) {
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1A);
//   Wire.write(0x05);
//   Wire.endTransmission(); 
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1B);
//   Wire.write(0x08);
//   Wire.endTransmission();
//   Wire.beginTransmission(0x68);
//   Wire.write(0x43);
//   Wire.endTransmission(); 
//   Wire.requestFrom(0x68,6);
//   int16_t GyroX=Wire.read()<<8 | Wire.read();
//   int16_t GyroY=Wire.read()<<8 | Wire.read();
//   int16_t GyroZ=Wire.read()<<8 | Wire.read();
//   RateRoll=(float)GyroX/65.5;
//   RatePitch=(float)GyroY/65.5;
//   RateYaw=(float)GyroZ/65.5;
// }
// void setup() {
//   Serial.begin(57600);
//   pinMode(13, OUTPUT);
//   digitalWrite(13, HIGH);
//   Wire.setClock(400000);
//   Wire.begin();
//   delay(250);
//   Wire.beginTransmission(0x68); 
//   Wire.write(0x6B);
//   Wire.write(0x00);
//   Wire.endTransmission();
//   for (RateCalibrationNumber=0;
//          RateCalibrationNumber<2000; 
//          RateCalibrationNumber ++) {
//     gyro_signals();
//     RateCalibrationRoll+=RateRoll;
//     RateCalibrationPitch+=RatePitch;
//     RateCalibrationYaw+=RateYaw;
//     delay(1);
//   }
//   RateCalibrationRoll/=2000;
//   RateCalibrationPitch/=2000;
//   RateCalibrationYaw/=2000;   
// }
// void loop() {
//   gyro_signals();
//   RateRoll-=RateCalibrationRoll;
//   RatePitch-=RateCalibrationPitch;
//   RateYaw-=RateCalibrationYaw;

//   Serial.print("[w_bx, w_by, w_bz] [°/s] = [");
//   Serial.print(RateCalibrationRoll);
//   Serial.print(", ");
//   Serial.print(RateCalibrationPitch);
//   Serial.print(", ");
//   Serial.print(RateCalibrationYaw);
//   Serial.println("]");

//   Serial.print("[wx, wy, wz] [°/s] = [");
//   Serial.print(RateRoll);
//   Serial.print(", ");
//   Serial.print(RatePitch);
//   Serial.print(", ");
//   Serial.print(RateYaw);
//   Serial.println("]");


//   delay(50);
// }




// // Part IV: MPU-6050 gyroscope
// /*
// The contents of this code and instructions are the intellectual property of Carbon Aeronautics. 
// The text and figures in this code and instructions are licensed under a Creative Commons Attribution - Noncommercial - ShareAlike 4.0 International Public Licence. 
// This license lets you remix, adapt, and build upon your work non-commercially, as long as you credit Carbon Aeronautics 
// (but not in any way that suggests that we endorse you or your use of the work) and license your new creations under the identical terms.
// This code and instruction is provided "As Is” without any further warranty. Neither Carbon Aeronautics or the author has any liability to any person or entity 
// with respect to any loss or damage caused or declared to be caused directly or indirectly by the instructions contained in this code or by 
// the software and hardware described in it. As Carbon Aeronautics has no control over the use, setup, assembly, modification or misuse of the hardware, 
// software and information described in this manual, no liability shall be assumed nor accepted for any resulting damage or injury. 
// By the act of copying, use, setup or assembly, the user accepts all resulting liability.

// 1.0  5 October 2022 -  initial release
// */

// #include <Arduino.h>
// #include <Wire.h>
// float RateRoll, RatePitch, RateYaw;
// void gyro_signals(void) {
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1A);
//   Wire.write(0x05);
//   Wire.endTransmission(); 
//   Wire.beginTransmission(0x68);
//   Wire.write(0x1B); 
//   Wire.write(0x8); 
//   Wire.endTransmission(); 
//   Wire.beginTransmission(0x68);
//   Wire.write(0x43);
//   Wire.endTransmission();
//   Wire.requestFrom(0x68,6);
//   int16_t GyroX=Wire.read()<<8 | Wire.read();
//   int16_t GyroY=Wire.read()<<8 | Wire.read();
//   int16_t GyroZ=Wire.read()<<8 | Wire.read();
//   RateRoll=(float)GyroX/65.5;
//   RatePitch=(float)GyroY/65.5;
//   RateYaw=(float)GyroZ/65.5;
// }
// void setup() {
//   Serial.begin(57600);
//   pinMode(13, OUTPUT);
//   digitalWrite(13, HIGH);
//   Wire.setClock(400000);
//   Wire.begin();
//   delay(250);
//   Wire.beginTransmission(0x68); 
//   Wire.write(0x6B);
//   Wire.write(0x00);
//   Wire.endTransmission();
// }
// void loop() {
//   gyro_signals();
  
//   Serial.print("[wx, wy, wz] [°/s] = [");
//   Serial.print(RateRoll);
//   Serial.print(", ");
//   Serial.print(RatePitch);
//   Serial.print(", ");
//   Serial.print(RateYaw);
//   Serial.println("]");

//   delay(50);
// }




