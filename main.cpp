#define _USE_MATH_DEFINES
#include <Novice.h>
#include <cmath>

#ifdef _DEBUG
#include <imgui.h>
#endif // _DEBUG

const char kWindowTitle[] = "LD2A_02_ワダ_ケイタ";

struct Vector3 {
	float x, y, z;
};

struct Matrix4x4 {
	float m[4][4];
};

struct Sphere {
	Vector3 center;
	float radius;
};

// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result;

	for (int column = 0; column < 4; column++) {
		for (int row = 0; row < 4; row++) {

			result.m[column][row] = (m1.m[column][0] * m2.m[0][row]) + (m1.m[column][1] * m2.m[1][row]) + (m1.m[column][2] * m2.m[2][row]) + (m1.m[column][3] * m2.m[3][row]);

		}
	}

	return result;
};

// 拡大縮小行列
Matrix4x4 MakeScaleMatrix(const Vector3& scale) {
	Matrix4x4 result = { 0 };

	result.m[0][0] = scale.x;
	result.m[1][1] = scale.y;
	result.m[2][2] = scale.z;
	result.m[3][3] = 1;

	return result;
};

// X軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 result = {
		{
			{ 1, 0, 0, 0 },
			{ 0, std::cos(radian), std::sin(radian), 0 },
			{ 0, -std::sin(radian), std::cos(radian), 0 },
			{ 0, 0, 0, 1 }
		},
	};

	return result;
};

// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result = {
		{
			{ std::cos(radian), 0, -std::sin(radian), 0 },
			{ 0, 1, 0, 0 },
			{ std::sin(radian), 0, std::cos(radian), 0 },
			{ 0, 0, 0, 1 }
		},
	};

	return result;
};

// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result = {
		{
			{ std::cos(radian), std::sin(radian), 0, 0 },
			{ -std::sin(radian), std::cos(radian), 0, 0 },
			{ 0, 0, 1, 0 },
			{ 0, 0, 0, 1 }
		},
	};

	return result;
};

// 平行移動行列
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
	Matrix4x4 result = { 0 };

	result.m[0][0] = 1;
	result.m[1][1] = 1;
	result.m[2][2] = 1;
	result.m[3][3] = 1;

	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;

	return result;
};

// 3次元アフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 result = { 0 };

	result = Multiply(MakeScaleMatrix(scale), Multiply(MakeRotateXMatrix(rotate.x), Multiply(MakeRotateYMatrix(rotate.y), MakeRotateZMatrix(rotate.z))));
	result = Multiply(result, MakeTranslateMatrix(translate));

	return result;
};

// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {
			{
				{ 1.0f / (aspectRatio * std::tan(fovY / 2.0f)), 0, 0, 0},
				{ 0, 1.0f / std::tan(fovY / 2.0f), 0, 0},
				{ 0, 0, farClip / (farClip - nearClip), 1},
				{ 0, 0, (-nearClip * farClip) / (farClip - nearClip), 0},
			},
	};
	return result;
};

// 正射影行列
Matrix4x4 MakeOrthographMatrix(float left, float top, float right, float bottom, float nearClip, float farClip) {
	Matrix4x4 result = {
		{
			{ 2.0f / (right - left), 0, 0, 0 },
			{ 0, 2.0f / (top - bottom), 0, 0 },
			{ 0, 0, 1.0f / (farClip - nearClip), 0 },
			{ (left + right) / (left - right), (top + bottom) / (bottom - top), nearClip / (nearClip - farClip), 1 },
		},
	};
	return result;
};

// ビューポート行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 result = {
				{
					{ width / 2.0f, 0, 0, 0},
					{ 0, -(height / 2.0f), 0, 0},
					{ 0, 0, maxDepth - minDepth, 0},
					{ left + (width / 2.0f), top + (height / 2.0f), minDepth, 1},
				},
	};
	return result;
};

// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result{};

	result.x = (vector.x * matrix.m[0][0]) + (vector.y * matrix.m[1][0]) + (vector.z * matrix.m[2][0]) + (1.0f * matrix.m[3][0]);
	result.y = (vector.x * matrix.m[0][1]) + (vector.y * matrix.m[1][1]) + (vector.z * matrix.m[2][1]) + (1.0f * matrix.m[3][1]);
	result.z = (vector.x * matrix.m[0][2]) + (vector.y * matrix.m[1][2]) + (vector.z * matrix.m[2][2]) + (1.0f * matrix.m[3][2]);
	float w = (vector.x * matrix.m[0][3]) + (vector.y * matrix.m[1][3]) + (vector.z * matrix.m[2][3]) + (1.0f * matrix.m[3][3]);
	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
};

// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m) {
	Matrix4x4 result;

	float x = 1 / (m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3]
		+ m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1]
		+ m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]

		- m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1]
		- m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3]
		- m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]

		- m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3]
		- m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1]
		- m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]

		+ m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1]
		+ m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3]
		+ m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]

		+ m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3]
		+ m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1]
		+ m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]

		- m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1]
		- m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3]
		- m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]

		- m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0]
		- m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0]
		- m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]

		+ m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0]
		+ m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0]
		+ m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0]
		);

	result.m[0][0] = ((m.m[1][1] * m.m[2][2] * m.m[3][3]) + (m.m[1][2] * m.m[2][3] * m.m[3][1]) + (m.m[1][3] * m.m[2][1] * m.m[3][2])
		- (m.m[1][3] * m.m[2][2] * m.m[3][1]) - (m.m[1][2] * m.m[2][1] * m.m[3][3]) - (m.m[1][1] * m.m[2][3] * m.m[3][2])) * x;

	result.m[0][1] = ((-m.m[0][1] * m.m[2][2] * m.m[3][3]) - (m.m[0][2] * m.m[2][3] * m.m[3][1]) - (m.m[0][3] * m.m[2][1] * m.m[3][2])
		+ (m.m[0][3] * m.m[2][2] * m.m[3][1]) + (m.m[0][2] * m.m[2][1] * m.m[3][3]) + (m.m[0][1] * m.m[2][3] * m.m[3][2])) * x;

	result.m[0][2] = ((m.m[0][1] * m.m[1][2] * m.m[3][3]) + (m.m[0][2] * m.m[1][3] * m.m[3][1]) + (m.m[0][3] * m.m[1][1] * m.m[3][2])
		- (m.m[0][3] * m.m[1][2] * m.m[3][1]) - (m.m[0][2] * m.m[1][1] * m.m[3][3]) - (m.m[0][1] * m.m[1][3] * m.m[3][2])) * x;

	result.m[0][3] = ((-m.m[0][1] * m.m[1][2] * m.m[2][3]) - (m.m[0][2] * m.m[1][3] * m.m[2][1]) - (m.m[0][3] * m.m[1][1] * m.m[2][2])
		+ (m.m[0][3] * m.m[1][2] * m.m[2][1]) + (m.m[0][2] * m.m[1][1] * m.m[2][3]) + (m.m[0][1] * m.m[1][3] * m.m[2][2])) * x;


	result.m[1][0] = ((-m.m[1][0] * m.m[2][2] * m.m[3][3]) - (m.m[1][2] * m.m[2][3] * m.m[3][0]) - (m.m[1][3] * m.m[2][0] * m.m[3][2])
		+ (m.m[1][3] * m.m[2][2] * m.m[3][0]) + (m.m[1][2] * m.m[2][0] * m.m[3][3]) + (m.m[1][0] * m.m[2][3] * m.m[3][2])) * x;

	result.m[1][1] = ((m.m[0][0] * m.m[2][2] * m.m[3][3]) + (m.m[0][2] * m.m[2][3] * m.m[3][0]) + (m.m[0][3] * m.m[2][0] * m.m[3][2])
		- (m.m[0][3] * m.m[2][2] * m.m[3][0]) - (m.m[0][2] * m.m[2][0] * m.m[3][3]) - (m.m[0][0] * m.m[2][3] * m.m[3][2])) * x;

	result.m[1][2] = ((-m.m[0][0] * m.m[1][2] * m.m[3][3]) - (m.m[0][2] * m.m[1][3] * m.m[3][0]) - (m.m[0][3] * m.m[1][0] * m.m[3][2])
		+ (m.m[0][3] * m.m[1][2] * m.m[3][0]) + (m.m[0][2] * m.m[1][0] * m.m[3][3]) + (m.m[0][0] * m.m[1][3] * m.m[3][2])) * x;

	result.m[1][3] = ((m.m[0][0] * m.m[1][2] * m.m[2][3]) + (m.m[0][2] * m.m[1][3] * m.m[2][0]) + (m.m[0][3] * m.m[1][0] * m.m[2][2])
		- (m.m[0][3] * m.m[1][2] * m.m[2][0]) - (m.m[0][2] * m.m[1][0] * m.m[2][3]) - (m.m[0][0] * m.m[1][3] * m.m[2][2])) * x;


	result.m[2][0] = ((m.m[1][0] * m.m[2][1] * m.m[3][3]) + (m.m[1][1] * m.m[2][3] * m.m[3][0]) + (m.m[1][3] * m.m[2][0] * m.m[3][1])
		- (m.m[1][3] * m.m[2][1] * m.m[3][0]) - (m.m[1][1] * m.m[2][0] * m.m[3][3]) - (m.m[1][0] * m.m[2][3] * m.m[3][1])) * x;

	result.m[2][1] = ((-m.m[0][0] * m.m[2][1] * m.m[3][3]) - (m.m[0][1] * m.m[2][3] * m.m[3][0]) - (m.m[0][3] * m.m[2][0] * m.m[3][1])
		+ (m.m[0][3] * m.m[2][1] * m.m[3][0]) + (m.m[0][1] * m.m[2][0] * m.m[3][3]) + (m.m[0][0] * m.m[2][3] * m.m[3][1])) * x;

	result.m[2][2] = ((m.m[0][0] * m.m[1][1] * m.m[3][3]) + (m.m[0][1] * m.m[1][3] * m.m[3][0]) + (m.m[0][3] * m.m[1][0] * m.m[3][1])
		- (m.m[0][3] * m.m[1][1] * m.m[3][0]) - (m.m[0][1] * m.m[1][0] * m.m[3][3]) - (m.m[0][0] * m.m[1][3] * m.m[3][1])) * x;

	result.m[2][3] = ((-m.m[0][0] * m.m[1][1] * m.m[2][3]) - (m.m[0][1] * m.m[1][3] * m.m[2][0]) - (m.m[0][3] * m.m[1][0] * m.m[2][1])
		+ (m.m[0][3] * m.m[1][1] * m.m[2][0]) + (m.m[0][1] * m.m[1][0] * m.m[2][3]) + (m.m[0][0] * m.m[1][3] * m.m[2][1])) * x;


	result.m[3][0] = ((-m.m[1][0] * m.m[2][1] * m.m[3][2]) - (m.m[1][1] * m.m[2][2] * m.m[3][0]) - (m.m[1][2] * m.m[2][0] * m.m[3][1])
		+ (m.m[1][2] * m.m[2][1] * m.m[3][0]) + (m.m[1][1] * m.m[2][0] * m.m[3][2]) + (m.m[1][0] * m.m[2][2] * m.m[3][1])) * x;

	result.m[3][1] = ((m.m[0][0] * m.m[2][1] * m.m[3][2]) + (m.m[0][1] * m.m[2][2] * m.m[3][0]) + (m.m[0][2] * m.m[2][0] * m.m[3][1])
		- (m.m[0][2] * m.m[2][1] * m.m[3][0]) - (m.m[0][1] * m.m[2][0] * m.m[3][2]) - (m.m[0][0] * m.m[2][2] * m.m[3][1])) * x;

	result.m[3][2] = ((-m.m[0][0] * m.m[1][1] * m.m[3][2]) - (m.m[0][1] * m.m[1][2] * m.m[3][0]) - (m.m[0][2] * m.m[1][0] * m.m[3][1])
		+ (m.m[0][2] * m.m[1][1] * m.m[3][0]) + (m.m[0][1] * m.m[1][0] * m.m[3][2]) + (m.m[0][0] * m.m[1][2] * m.m[3][1])) * x;

	result.m[3][3] = ((m.m[0][0] * m.m[1][1] * m.m[2][2]) + (m.m[0][1] * m.m[1][2] * m.m[2][0]) + (m.m[0][2] * m.m[1][0] * m.m[2][1])
		- (m.m[0][2] * m.m[1][1] * m.m[2][0]) - (m.m[0][1] * m.m[1][0] * m.m[2][2]) - (m.m[0][0] * m.m[1][2] * m.m[2][1])) * x;

	return result;
};


void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSobdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSobdivision);
	Vector3 linePos[2]{};
	// 奥から手前への線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSobdivision; ++xIndex) {
		linePos[0] = { -kGridHalfWidth, 0.0f, kGridHalfWidth - (kGridEvery * xIndex) };
		linePos[1] = { kGridHalfWidth, 0.0f, kGridHalfWidth - (kGridEvery * xIndex) };

		for (uint32_t i = 0; i < 2; i++) {
			Vector3 ndcVertex = Transform(linePos[i], viewProjectionMatrix);
			linePos[i] = Transform(ndcVertex, viewportMatrix);
		}
		if (xIndex == 5) {
			Novice::DrawLine(int(linePos[0].x), int(linePos[0].y), int(linePos[1].x), int(linePos[1].y), 0x000000FF);
		}
		else {
			Novice::DrawLine(int(linePos[0].x), int(linePos[0].y), int(linePos[1].x), int(linePos[1].y), 0xAAAAAAFF);
		}
	}
	// 左から右への線を順々に引いていく
	for (uint32_t zIndex = 0; zIndex <= kSobdivision; ++zIndex) {
		linePos[0] = { kGridHalfWidth - (kGridEvery * zIndex), 0.0f, -kGridHalfWidth };
		linePos[1] = { kGridHalfWidth - (kGridEvery * zIndex), 0.0f, kGridHalfWidth };

		for (uint32_t i = 0; i < 2; i++) {
			Vector3 ndcVertex = Transform(linePos[i], viewProjectionMatrix);
			linePos[i] = Transform(ndcVertex, viewportMatrix);
		}
		if (zIndex == 5) {
			Novice::DrawLine(int(linePos[0].x), int(linePos[0].y), int(linePos[1].x), int(linePos[1].y), 0x000000FF);
		}
		else {
			Novice::DrawLine(int(linePos[0].x), int(linePos[0].y), int(linePos[1].x), int(linePos[1].y), 0xAAAAAAFF);
		}
	}
}

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 10; // 分割数
	const float kLonEvery = (2.0f * (float)M_PI) / kSubdivision; // 経度分割1つ分の角度
	const float kLatEvery = (float)M_PI / kSubdivision; // 緯度分割1つ分の角度

	// 緯度の方向に分割 -π/2 ～ π/2
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -(float)M_PI / 2.0f + kLatEvery * latIndex; // 現在の緯度

		// 緯度方向に分割 0 ～ 2π
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;
			// world座標系でのa,b,cを求める
			Vector3 a, b, c;

			a = { sphere.radius * cosf(lat) * cosf(lon), sphere.radius * sinf(lat), sphere.radius * cosf(lat) * sin(lon) };
			b = { sphere.radius * cosf(lat + kLatEvery) * cosf(lon), sphere.radius * sinf(lat + kLatEvery), sphere.radius * cosf(lat + kLatEvery) * sin(lon) };
			c = { sphere.radius * cosf(lat) * cosf(lon + kLonEvery), sphere.radius * sinf(lat), sphere.radius * cosf(lat) * sin(lon + kLonEvery) };

			// a,b,cをScreen座標系まで変換
			Vector3 ndcVertex = Transform(a, viewProjectionMatrix);
			a = Transform(ndcVertex, viewportMatrix);

			ndcVertex = Transform(b, viewProjectionMatrix);
			b = Transform(ndcVertex, viewportMatrix);

			ndcVertex = Transform(c, viewProjectionMatrix);
			c = Transform(ndcVertex, viewportMatrix);

			// ab,bcで線を引く
			Novice::DrawLine((int)a.x, (int)a.y, (int)b.x, (int)b.y, color);
			Novice::DrawLine((int)a.x, (int)a.y, (int)c.x, (int)c.y, color);
		}
	}
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	const int kWindowWidth = 1280;
	const int kWindowHeight = 720;

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	Vector3 scale{ 1.0f,1.0f,1.0f };
	Vector3 rotate{};
	Vector3 translate{};

	Vector3 sphereScale{ 1.0f,1.0f,1.0f };
	Vector3 sphereTranslate{};
	Vector3 sphereRotate{};

	Vector3 cameraScale{ 1.0f,1.0f,1.0f };
	Vector3 cameraPosition{ 0.0f,1.9f,-6.49f };
	Vector3 cameraRotation{ 0.26f,0.0f,0.0f };

	Sphere sphere{ {0,0,0}, {1.0f} };


	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///
		Matrix4x4 sphereWorldMatrix = MakeAffineMatrix(sphereScale, sphereRotate, sphereTranslate);
		Matrix4x4 worldMatrix = MakeAffineMatrix(scale, rotate, translate);

		Matrix4x4 cameraMatrix = MakeAffineMatrix(cameraScale, cameraRotation, cameraPosition);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);

		Matrix4x4 sphereWorldViewProjectionMatrix = Multiply(sphereWorldMatrix, Multiply(viewMatrix, projectionMatrix));

		Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));

		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		ImGui::Begin("camera");
		ImGui::DragFloat3("scale", &cameraScale.x, 0.01f);
		ImGui::DragFloat3("position", &cameraPosition.x, 0.01f);
		ImGui::DragFloat3("rotate", &cameraRotation.x, 0.01f);
		ImGui::End();

		ImGui::Begin("sphere");
		ImGui::DragFloat3("scale", &sphereScale.x, 0.01f);
		ImGui::DragFloat3("rotate", &sphereRotate.x, 0.01f);
		ImGui::DragFloat3("translate", &sphereTranslate.x, 0.01f);
		ImGui::DragFloat("radius", &sphere.radius, 0.01f);
		ImGui::End();

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		DrawGrid(worldViewProjectionMatrix, viewportMatrix);
		DrawSphere(sphere, sphereWorldViewProjectionMatrix, viewportMatrix, 0xFFFFFFFF);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
