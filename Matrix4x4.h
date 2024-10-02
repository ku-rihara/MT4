#pragma once
#include"Vector3.h"

class Matrix4x4 {

public:
	float m[4][4];
public:

	Matrix4x4();
	Matrix4x4(float a1,float a2,float a3, float a4,float b1, float b2, float b3,float b4, float c1, float c2, float c3,float c4, float d1, float d2, float d3, float d4);
	/*Matrix4x4();*/

	Matrix4x4 operator+(const Matrix4x4& obj);//加算

	Matrix4x4 operator-(const Matrix4x4& obj);//減算

	Matrix4x4 operator*(const Matrix4x4& obj);//積

	static	Matrix4x4 Transpose(const Matrix4x4& m);//転置行列

	static	Matrix4x4 MakeIdentity4x4();//単位行列

	Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2);

	static	Matrix4x4 MakeTranslateMatrix(const Vector3& translate);//平行移動行列

	static	Matrix4x4 MakeScaleMatrix(const Vector3& scale);//拡大縮小行列

	static	Matrix4x4 MakeRotateXMatrix(float radian);// X軸回転行列

	static	Matrix4x4 MakeRotateYMatrix(float radian);// Y軸回転行列

	static	Matrix4x4 MakeRotateZMatrix(float radian);// Z軸回転行列

	static	Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate);//アフィン変換

	static	Matrix4x4 Inverse(const Matrix4x4& m);//逆行列

	static	Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip);//透視投影行列

	static	Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip);//正射影行列

	static	Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth);//ビューポート行列

	static	Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix);//スクリーン変換
	static Vector3 ScreenTransform(const Vector3& local, const Matrix4x4& worldviewprojection, const Matrix4x4& viewport);

	static	void MatrixScreenPrintf(int x, int y, Matrix4x4& matrix, const char* label);
};
