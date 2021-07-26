#include <cstdlib> //pour faire une pause sur l'affichage
#include <iostream>
#include <cmath>

#include <FreeImage.h>

using namespace std;

class Transformation {
public:

	Transformation() {
		mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;//premiere ligne
		mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = 0;//deuxieme ligne
		mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1; mat[2][3] = 0;//troisieme ligne
		mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;//quatrieme ligne
	}

	virtual void print() {
		cout << "class generique : Identite" << endl;
	}

	Transformation* compose(Transformation Tr) { //compose_gauche

		float tmp[4][4];

		for (int i = 0; i < 4; i++)

			for (int j = 0; j < 4; j++)
			{
				tmp[i][j] = 0;
				for (int k = 0; k < 4; k++)
					tmp[i][j] += Tr.mat[i][k] * mat[k][j];
			}

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				mat[i][j] = tmp[i][j];

		return this;
	}

	float mat[4][4];
};


class Point {
public:
	Point() {
		x = 0;
		y = 0;
		z = 0;
		w = 1;
	}

	Point(float x, float y, float z) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = 1;
	}
	Point(float* V) {
		x = V[0];
		y = V[1];
		z = V[2];
		w = 1;

	}
	Point* Transform(Transformation& trsf) {

		float tx = trsf.mat[0][0] * x + trsf.mat[0][1] * y + trsf.mat[0][2] * z + trsf.mat[0][3] * w;
		float ty = trsf.mat[1][0] * x + trsf.mat[1][1] * y + trsf.mat[1][2] * z + trsf.mat[1][3] * w;
		float tz = trsf.mat[2][0] * x + trsf.mat[2][1] * y + trsf.mat[2][2] * z + trsf.mat[2][3] * w;
		float tw = trsf.mat[3][0] * x + trsf.mat[3][1] * y + trsf.mat[3][2] * z + trsf.mat[3][3] * w;


		this->x = tx / tw;
		this->y = ty / tw;
		this->z = tz / tw;
		this->w = 1;

		//cout << "tw  " << tw << endl;

		//cout << "x " << x << "  y  " << y << "  z  " << z << endl;

		return this;
	}

	float x;
	float y;
	float z;
	float w;

	void print() {

		//cout << "x= " << x << " , y= " << y << " , z= " << z << endl;
	}

};

class Vector {
public:
	Vector() {
		x = 0;
		y = 0;
		z = 0;
		w = 0;
	}

	Vector(float x, float y, float z) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = 0;
	}

	float x;
	float y;
	float z;
	float w;

	void print() {

		//cout << "x= " << x << " , y= " << y << " , z= " << z << endl;
	}

	float norm() {
		return sqrt(x * x + y * y + z * z);
	}

	Vector normalize() {
		float norme = norm();
		return Vector(x / norme, y / norme, z / norme);
	}

	Vector operator- () {
		return Vector(-x, -y, -z);
	}

	Vector operator^ (Vector V) {

		double v1x = this->x, v1y = this->y, v1z = this->z;
		double v2x = V.x, v2y = V.y, v2z = V.z;

		return Vector((v1y * v2z) - (v1z * v2y),
			(v1z * v2x) - (v1x * v2z),
			(v1x * v2y) - (v1y * v2x));

	}
};

class Rotation_Z : public Transformation {
public:

	Rotation_Z(float ang) {

		angle = ang / 180 * 3.14159;

		mat[0][0] = cosf(angle); mat[0][1] = -sinf(angle);  mat[0][2] = 0; mat[0][3] = 0;//premiere ligne
		mat[1][0] = sinf(angle); mat[1][1] = cosf(angle);  mat[1][2] = 0; mat[1][3] = 0;//deuxieme ligne
		mat[2][0] = 0;           mat[2][1] = 0;            mat[2][2] = 1; mat[2][3] = 0;//troisieme ligne
		mat[3][0] = 0;           mat[3][1] = 0;            mat[3][2] = 0; mat[3][3] = 1;//quatrieme ligne
	}
	float angle;

	virtual void print() {
		cout << " je suis une rotation, d'un angle " << angle << endl;
	}

};


class Rotation_X : public Transformation {
public:

	Rotation_X(float ang) {

		angle = ang / 180 * 3.14159;

		mat[0][0] = 1;  mat[0][1] = 0;  mat[0][2] = 0;            mat[0][3] = 0;//premiere ligne
		mat[1][0] = 0;  mat[1][1] = cosf(angle);  mat[1][2] = -sinf(angle);  mat[1][3] = 0;//deuxieme ligne
		mat[2][0] = 0;  mat[2][1] = sinf(angle);  mat[2][2] = cosf(angle);  mat[2][3] = 0;//troisieme ligne
		mat[3][0] = 0;  mat[3][1] = 0;            mat[3][2] = 0;            mat[3][3] = 1;//quatrieme ligne
	}
	float angle;

	virtual void print() {
		cout << " je suis une rotation, d'un angle " << angle << endl;
	}

};

class Rotation_Y : public Transformation {
public:

	Rotation_Y(float ang) {

		angle = ang / 180 * 3.14159;

		mat[0][0] = cosf(angle);  mat[0][1] = 0;  mat[0][2] = sinf(angle);  mat[0][3] = 0;//premiere ligne
		mat[1][0] = 0;            mat[1][1] = 1;  mat[1][2] = 0;            mat[1][3] = 0;//deuxieme ligne
		mat[2][0] = -sinf(angle); mat[2][1] = 0;  mat[2][2] = cosf(angle); mat[2][3] = 0;//troisieme ligne
		mat[3][0] = 0;            mat[3][1] = 0;  mat[3][2] = 0;            mat[3][3] = 1;//quatrieme ligne
	}
	float angle;

	virtual void print() {
		cout << " je suis une rotation, d'un angle " << angle << endl;
	}

};

class Translation : public Transformation {
public:

	Translation(Vector V) {

		mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = V.x; //premiere ligne
		mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = V.y; //deuxieme ligne
		mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 0; mat[2][3] = V.z;   //troisieme ligne
		mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;   //troisieme ligne

		this->V = V;

	}

	virtual void print() {
		cout << " je suis une translation, d'un vecteur "; V.print();
	}

	Vector V;

};

class Viewport : public Transformation {
public:

	Viewport(int n_x, int n_y) {

		mat[0][0] = float(n_x) / 2.0; mat[0][1] = 0;              mat[0][2] = 0; mat[0][3] = float(n_x - 1) / 2.0; //premiere ligne
		mat[1][0] = 0;               mat[1][1] = float(n_y) / 2.0; mat[1][2] = 0; mat[1][3] = float(n_y - 1) / 2.0; //deuxieme ligne
		mat[2][0] = 0;               mat[2][1] = 0;              mat[2][2] = 1; mat[2][3] = 0;   //troisieme ligne
		mat[3][0] = 0;               mat[3][1] = 0;              mat[3][2] = 0; mat[3][3] = 1;   //troisieme ligne

		res_x = n_x;
		res_y = n_y;
	}

	virtual void print() {
		cout << " je suis Viewport res_x= " << res_x << "res_y= " << res_y << endl;
	}

	int res_x;
	int res_y;
};

class Perspective : public Transformation {
public:

	Perspective(float ang, float ratio, float n, float f) {

		float angle = ang / 180 * 3.14159;

		mat[0][0] = ratio / tanf(angle / 2.0); mat[0][1] = 0;                  mat[0][2] = 0;           mat[0][3] = 0; //premiere ligne
		mat[1][0] = 0;                  mat[1][1] = 1 / tanf(angle / 2.0);  mat[1][2] = 0;           mat[1][3] = 0; //deuxieme ligne
		mat[2][0] = 0;                  mat[2][1] = 0;                      mat[2][2] = (f + n) / (n - f); mat[2][3] = 2 * f * n / (f - n);   //troisieme ligne
		mat[3][0] = 0;                  mat[3][1] = 0;                      mat[3][2] = 1;           mat[3][3] = 0;   //troisieme ligne

		this->angle = angle;
		this->ratio = ratio;
		this->n = n;
		this->f = f;
	}

	virtual void print() {
		cout << " je suis Perspective : angle= " << angle << " ratio= " << ratio << "n= " << n << "f= " << f << endl;
	}

	float angle; float ratio; float n; float f;
};

class LookAt : public Transformation {
public:

	LookAt(Point eye, Vector dir, Vector up) {

		this->eye = eye;
		this->dir = dir;
		this->up = up;

		z_cam = -dir.normalize();

		x_cam = (up ^ z_cam).normalize();

		y_cam = z_cam ^ x_cam;

		float mat1[4][4];
		float mat2[4][4];


		mat1[0][0] = x_cam.x; mat1[0][1] = x_cam.y; mat1[0][2] = x_cam.z; mat1[0][3] = 0;
		mat1[1][0] = y_cam.x; mat1[1][1] = y_cam.y; mat1[1][2] = y_cam.z; mat1[1][3] = 0;
		mat1[2][0] = z_cam.x; mat1[2][1] = z_cam.y; mat1[2][2] = z_cam.z; mat1[2][3] = 0;
		mat1[3][0] = 0;       mat1[3][1] = 0;       mat1[3][2] = 0;       mat1[3][3] = 1;

		mat2[0][0] = 1; mat2[0][1] = 0; mat2[0][2] = 0; mat2[0][3] = -eye.x;
		mat2[1][0] = 0; mat2[1][1] = 1; mat2[1][2] = 0; mat2[1][3] = -eye.y;
		mat2[2][0] = 0; mat2[2][1] = 0; mat2[2][2] = 1; mat2[2][3] = -eye.z;
		mat2[3][0] = 0; mat2[3][1] = 0; mat2[3][2] = 0; mat2[3][3] = 1;

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				mat[i][j] = 0;
				for (int k = 0; k < 4; k++) {
					mat[i][j] = mat[i][j] + mat1[i][k] * mat2[k][j];
				}
			}

	}

	virtual void print() {
		cout << "Je suis LookAt : " << endl;
		cout << "eye ";  eye.print();
		cout << "dir "; dir.print();
		cout << "up "; up.print();
		cout << "x_cam "; x_cam.print();
		cout << "y_cam "; y_cam.print();
		cout << "z_cam "; z_cam.print();
	}

	Point eye; Vector dir; Vector up; Vector x_cam; Vector  y_cam; Vector z_cam;
};


//TODO: A optimiser ...
int DDA(Point p1, Point p2, FIBITMAP* bitmap, RGBQUAD color) {

	Point start_p;
	Point end_p;

	if (p1.x <= p2.x)
	{
		start_p.x = p1.x;
		start_p.y = p1.y;

		end_p.x = p2.x;
		end_p.y = p2.y;

	}
	else {
		start_p.x = p2.x;
		start_p.y = p2.y;

		end_p.x = p1.x;
		end_p.y = p1.y;
	}

	//cout << "start ";  start_p.print();
	//cout << "end "; end_p.print();

	int delta_x = end_p.x - start_p.x;
	int delta_y = end_p.y - start_p.y;

	float m = (double)delta_y / (double)delta_x;


	FreeImage_SetPixelColor(bitmap, floor(start_p.x + 0.5), floor(start_p.y + 0.5), &color);

	if (delta_x == 0) {

		if (start_p.y < end_p.y)
			for (int y = start_p.y + 1; y < end_p.y; y++)
				FreeImage_SetPixelColor(bitmap, start_p.x, y, &color);
		else
			for (int y = end_p.y + 1; y < start_p.y; y++)
				FreeImage_SetPixelColor(bitmap, start_p.x, y, &color);

	}
	else if (delta_y == 0)
	{
		if (start_p.x < end_p.x)
			for (int x = start_p.x + 1; x < end_p.x; x++)
				FreeImage_SetPixelColor(bitmap, x, start_p.y, &color);
		else
			for (int x = end_p.x + 1; x < start_p.x; x++)
				FreeImage_SetPixelColor(bitmap, x, start_p.y, &color);


	}
	else if (abs(m) < 1) {

		//cout << "m < 1 " << endl;
		//cout << "m =  " << m << endl;

		float y = start_p.y;

		for (int x = start_p.x + 1; x < end_p.x; x++)
		{
			y = y + m;
			FreeImage_SetPixelColor(bitmap, x, floor(y + 0.5), &color);

		}
	}
	else {
		//cout << "m =  " << m << endl;
		float x = start_p.x;

		if (start_p.y < end_p.y)
			for (int y = start_p.y + 1; y < end_p.y; y++)
			{
				x = x + 1.0 / m;
				FreeImage_SetPixelColor(bitmap, floor(x + 0.5), y, &color);

			}
		else
			for (int y = start_p.y - 1; y > end_p.y; y--)
			{
				x = x - 1.0 / m;
				FreeImage_SetPixelColor(bitmap, floor(x + 0.5), y, &color);

			}

	}

	FreeImage_SetPixelColor(bitmap, floor(end_p.x + 0.5), floor(end_p.y + 0.5), &color);

	return 0;
}

void Dessiner_Quad(Point p1, Point p2, Point p3, Point p4, FIBITMAP * bitmap, RGBQUAD color) {
	DDA(p1, p2, bitmap, color);
	DDA(p2, p3, bitmap, color);
	DDA(p3, p4, bitmap, color);
	DDA(p4, p1, bitmap, color);
}

int main() {

	FreeImage_Initialise();
	cout << "FreeImage " << FreeImage_GetVersion() << endl << endl;
	//FreeImage_DeInitialise();

	int largeur = 2000;
	int hauteur = 2000;

	int BPP = 24; //Bits par pixel (8 bits pour chaque caanal RGB

	FIBITMAP* bitmap = FreeImage_Allocate(largeur, hauteur, BPP);

	RGBQUAD color;
	color.rgbRed = 255.0;
	color.rgbGreen = 255.0;
	color.rgbBlue = 255.0;

	LookAt cam(Point(0, 13, 0), Vector(0, -1, 0), Vector(0, 1, 1));
	Perspective persp(45, 1, -0.1, -100);
	Viewport viewport(largeur, hauteur);

	//Sera appliqu�e � tous les sommets
	Transformation Tf;
	Tf.compose(Rotation_Z(10))->compose(cam)->compose(persp)->compose(viewport);

	/*//cube
	float v[8][3];
	v[0][0] = v[1][0] = v[2][0] = v[3][0] = -1;
	v[4][0] = v[5][0] = v[6][0] = v[7][0] = 1;
	v[0][1] = v[1][1] = v[4][1] = v[5][1] = -1;
	v[2][1] = v[3][1] = v[6][1] = v[7][1] = 1;
	v[0][2] = v[3][2] = v[4][2] = v[7][2] = 1;
	v[1][2] = v[2][2] = v[5][2] = v[6][2] = -1;

	int faces[6][4] = {  /* les indices des sommets des 6 faces du cube.
		 { 0, 1, 2, 3 },{ 3, 2, 6, 7 },{ 7, 6, 5, 4 },{ 4, 5, 1, 0 },{ 5, 6, 2, 1 },{ 7, 4, 0, 3 } };

	for (int i = 0; i < 6; i++) {

		Point p1(&v[faces[i][0]][0]), p2(&v[faces[i][1]][0]), p3(&v[faces[i][2]][0]), p4(&v[faces[i][3]][0]);

		p1.Transform(Tf);
		p2.Transform(Tf);
		p3.Transform(Tf);
		p4.Transform(Tf);

		Dessiner_Quad(p1, p2, p3, p4, bitmap, color);

	}   */
	// fifinir les pointde la sphere ***************************************************
	double x, y, z;
	Point mat[72][36];
	float rayon = 4;
	for (int tetha = 0; tetha < 360; tetha = tetha + 5)
		for (int fi = 0; fi < 180; fi = fi + 5)
		{
			x = rayon * sin(fi * 3.14 / 180) * cos(tetha * 3.14 / 180);
			y = rayon * sin(fi * 3.14 / 180) * sin(tetha * 3.14 / 180);
			z = rayon * cos(fi * 3.14 / 180);
			Point p(x, y, z);
			mat[tetha / 5][fi / 5] = p;
		}
	// la tranformation de point (transformation-> cam -> perspective -> view port)
	for (int i = 0; i < 72; i++)
		for (int j = 0; j < 36; j++)
		{
			mat[i][j].Transform(Tf);
		}
	//dessin de la sphere
	for (int i = 0; i < 72; i++)
		for (int j = 0; j < 35; j++)//36
		{
			DDA(mat[i][j], mat[i][j + 1], bitmap, color);
		}
	for (int i = 0; i < 36; i++)//72  j
		for (int j = 0; j < 71; j++)//72  i
		{
			DDA(mat[j][i], mat[j + 1][i], bitmap, color);
			//DDA(mat[i][j], mat[i][j + 1], bitmap, color);
		}
	for (int i = 0; i < 71; i++)//72  j
		for (int j = 0; j < 35; j++)//72  i
		{
			DDA(mat[i][j], mat[i+1][j+1], bitmap, color);
			//DDA(mat[i][j], mat[i][j + 1], bitmap, color);
		}


	if (FreeImage_Save(FIF_PNG, bitmap, "tp1.png", 0)) {
		cout << "Image successfully saved!" << endl;
		FreeImage_Unload(bitmap);
	}

	FreeImage_DeInitialise();

	system("PAUSE");

	return 0;

}