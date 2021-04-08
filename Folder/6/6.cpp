#include <iostream>
#include <ctime>
using namespace std;

int ch = 1;

//копируем в линейно размещенную матрицу a квадрант матрицы b
//начинаем с элемента b[ib][[jb]
//n - размерность матрицы b
void copyM(int* a, int* b, int ib, int jb, int n) {
	int i, j, k;
	int imax = ib + n / 2;
	int jmax = jb + n / 2;

	for (k = 0, i = ib; i < imax; i++) {
		for (j = jb; j < jmax; j++) {
			a[k++] = b[i * n + j];
		}
	}
}

void copy(int* a, int** b, int ib, int jb, int n) {
	int i, j, k;
	int imax = ib + n / 2;
	int jmax = jb + n / 2;

	for (k = 0, i = ib; i < imax; i++) {
		for (j = jb; j < jmax; j++) {
			a[k++] = b[i][j];
		}
	}
}

//копируем в квадрант матрицы a линейно размещенную матрицу b
//начинаем с элемента a[ia][[ja]
//n - размерность матрицы a
void copybackM(int* a, int ia, int ja, int* b, int n) {
	int i, j, k;

	int imax = ia + n / 2;
	int jmax = ja + n / 2;

	for (k = 0, i = ia; i < imax; i++) {
		for (j = ja; j < jmax; j++) {
			a[i * n + j] = b[k++];
		}
	}
}

void copyback(int** a, int ia, int ja, int* b, int n) {
	int i, j, k;

	int imax = ia + n / 2;
	int jmax = ja + n / 2;

	for (k = 0, i = ia; i < imax; i++) {
		for (j = ja; j < jmax; j++) {
			a[i][j] = b[k++];
		}
	}
}

//складываем линейно размещенные матрицы c = a + b
void add(int* c, int* a, int* b, int n) {
	for (int i = 0; i < n * n; i++)
		c[i] = a[i] + b[i];
}

//вычитаем линейно размещенные матрицы c = a - b
void sub(int* c, int* a, int* b, int n) {
	for (int i = 0; i < n * n; i++)
		c[i] = a[i] - b[i];
}

//обычное матричное умножение c = a * b
void mul_normalM(int* c, int* a, int* b, int n) {
	int i, j, k;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			c[i * n + j] = 0;
			for (k = 0; k < n; k++)
				c[i * n + j] += a[i * n + k] * b[k * n + j];
		}
	}
}

void mul_normal(int** c, int** a, int** b, int n) {
	int i, j, k;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			c[i][j] = 0;
			for (k = 0; k < n; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
}

//умножение алгоритмом Винограда - Штрассена (модификация алгоритма Штрассена)
//c = a * b
void mulM(int* c, int* a, int* b, int n) {
	if (n <= 2) {
		mul_normalM(c, a, b, n);
	} else {
		int h = n / 2;                    
		int* M = new int[h * h * 29];       						

		copyM(&M[0], a, 0, 0, n);                    	 //M[0] = A11
		copyM(&M[h * h], a, 0, h, n);                  	 //M[1] = A12
		copyM(&M[2 * h * h], a, h, 0, n);                //M[2] = A21
		copyM(&M[3 * h * h], a, h, h, n);                //M[3] = A22
		copyM(&M[4 * h * h], b, 0, 0, n);                //M[4] = B11
		copyM(&M[5 * h * h], b, 0, h, n);                //M[5] = B12
		copyM(&M[6 * h * h], b, h, 0, n);                //M[6] = B21
		copyM(&M[7 * h * h], b, h, h, n);                //M[7] = B22
		add(&M[8 * h * h], &M[2 * h * h], &M[3 * h * h], h);    //M[8] = S1 = A21 + A22
		sub(&M[9 * h * h], &M[8 * h * h], &M[0], h);        	//M[9] = S2 = S1 - A11
		sub(&M[10 * h * h], &M[0], &M[2 * h * h], h);       	//M[10] = S3 = A11 - A21
		sub(&M[11 * h * h], &M[h * h], &M[9 * h * h], h);     	//M[11] = S4 = A12 - S2
		sub(&M[12 * h * h], &M[5 * h * h], &M[4 * h * h], h);   //M[12] = S5 = B12 - B11
		sub(&M[13 * h * h], &M[7 * h * h], &M[12 * h * h], h);  //M[13] = S6 = B22 - S5
		sub(&M[14 * h * h], &M[7 * h * h], &M[5 * h * h], h);   //M[14] = S7 = B22 - B12
		sub(&M[15 * h * h], &M[13 * h * h], &M[6 * h * h], h);  //M[15] = S8 = S6 - B21
		
		mulM(&M[16 * h * h], &M[9 * h * h], &M[13 * h * h], h);  	//M[16] = P1 = S2 * S6
		mulM(&M[17 * h * h], &M[0], &M[4 * h * h], h);       		//M[17] = P2 = A11 * B11
		mulM(&M[18 * h * h], &M[h * h], &M[6 * h * h], h);     		//M[18] = P3 = A12 * B21
		mulM(&M[19 * h * h], &M[10 * h * h], &M[14 * h * h], h); 	//M[19] = P4 = S3 * S7
		mulM(&M[20 * h * h], &M[8 * h * h], &M[12 * h * h], h);  	//M[20] = P5 = S1 * S5
		mulM(&M[21 * h * h], &M[11 * h * h], &M[7 * h * h], h);  	//M[21] = P6 = S4 * B22
		mulM(&M[22 * h * h], &M[3 * h * h], &M[15 * h * h], h);  	//M[22] = P7 = A22 * S8

		add(&M[23 * h * h], &M[16 * h * h], &M[17 * h * h], h); 	//M[23] = T1 = P1 + P2
		add(&M[24 * h * h], &M[23 * h * h], &M[19 * h * h], h); 	//M[24] = T2 = T1 + P4
		add(&M[25 * h * h], &M[17 * h * h], &M[18 * h * h], h); 	//M[25] = C11 = P2 + P3
		add(&M[26 * h * h], &M[23 * h * h], &M[20 * h * h], h); 	//M[26] = C12 = T1 + P5
		add(&M[26 * h * h], &M[26 * h * h], &M[21 * h * h], h); 	//M[26] = C12 += P6
		sub(&M[27 * h * h], &M[24 * h * h], &M[22 * h * h], h); 	//M[27] = C21 = T2 - P7
		add(&M[28 * h * h], &M[24 * h * h], &M[20 * h * h], h); 	//M[28] = C22 = T2 + P5

	//копируем результат
		copybackM(c, 0, 0, &M[25 * h * h], n);           //C11 = M[25]
		copybackM(c, 0, h, &M[26 * h * h], n);           //C12 = M[26]
		copybackM(c, h, 0, &M[27 * h * h], n);           //C21 = M[27]
		copybackM(c, h, h, &M[28 * h * h], n);           //C22 = M[28]
		delete M;
	}
}

void mul(int** c, int** a, int** b, int n) {
	if (n <= 2) {
		mul_normal(c, a, b, n);
	} else {
		int h = n / 2;                    
		int* M = new int[h * h * 29];       						

		copyM(&M[0], a, 0, 0, n);                    	 //M[0] = A11
		copyM(&M[h * h], a, 0, h, n);                  	 //M[1] = A12
		copyM(&M[2 * h * h], a, h, 0, n);                //M[2] = A21
		copyM(&M[3 * h * h], a, h, h, n);                //M[3] = A22
		copyM(&M[4 * h * h], b, 0, 0, n);                //M[4] = B11
		copyM(&M[5 * h * h], b, 0, h, n);                //M[5] = B12
		copyM(&M[6 * h * h], b, h, 0, n);                //M[6] = B21
		copyM(&M[7 * h * h], b, h, h, n);                //M[7] = B22
		add(&M[8 * h * h], &M[2 * h * h], &M[3 * h * h], h);    //M[8] = S1 = A21 + A22
		sub(&M[9 * h * h], &M[8 * h * h], &M[0], h);        	//M[9] = S2 = S1 - A11
		sub(&M[10 * h * h], &M[0], &M[2 * h * h], h);       	//M[10] = S3 = A11 - A21
		sub(&M[11 * h * h], &M[h * h], &M[9 * h * h], h);     	//M[11] = S4 = A12 - S2
		sub(&M[12 * h * h], &M[5 * h * h], &M[4 * h * h], h);   //M[12] = S5 = B12 - B11
		sub(&M[13 * h * h], &M[7 * h * h], &M[12 * h * h], h);  //M[13] = S6 = B22 - S5
		sub(&M[14 * h * h], &M[7 * h * h], &M[5 * h * h], h);   //M[14] = S7 = B22 - B12
		sub(&M[15 * h * h], &M[13 * h * h], &M[6 * h * h], h);  //M[15] = S8 = S6 - B21
		
		mulM(&M[16 * h * h], &M[9 * h * h], &M[13 * h * h], h);  	//M[16] = P1 = S2 * S6
		mulM(&M[17 * h * h], &M[0], &M[4 * h * h], h);       		//M[17] = P2 = A11 * B11
		mulM(&M[18 * h * h], &M[h * h], &M[6 * h * h], h);     		//M[18] = P3 = A12 * B21
		mulM(&M[19 * h * h], &M[10 * h * h], &M[14 * h * h], h); 	//M[19] = P4 = S3 * S7
		mulM(&M[20 * h * h], &M[8 * h * h], &M[12 * h * h], h);  	//M[20] = P5 = S1 * S5
		mulM(&M[21 * h * h], &M[11 * h * h], &M[7 * h * h], h);  	//M[21] = P6 = S4 * B22
		mulM(&M[22 * h * h], &M[3 * h * h], &M[15 * h * h], h);  	//M[22] = P7 = A22 * S8

		add(&M[23 * h * h], &M[16 * h * h], &M[17 * h * h], h); 	//M[23] = T1 = P1 + P2
		add(&M[24 * h * h], &M[23 * h * h], &M[19 * h * h], h); 	//M[24] = T2 = T1 + P4
		add(&M[25 * h * h], &M[17 * h * h], &M[18 * h * h], h); 	//M[25] = C11 = P2 + P3
		add(&M[26 * h * h], &M[23 * h * h], &M[20 * h * h], h); 	//M[26] = C12 = T1 + P5
		add(&M[26 * h * h], &M[26 * h * h], &M[21 * h * h], h); 	//M[26] = C12 += P6
		sub(&M[27 * h * h], &M[24 * h * h], &M[22 * h * h], h); 	//M[27] = C21 = T2 - P7
		add(&M[28 * h * h], &M[24 * h * h], &M[20 * h * h], h); 	//M[28] = C22 = T2 + P5

	//копируем результат
		copybackM(c, 0, 0, &M[25 * h * h], n);           //C11 = M[25]
		copybackM(c, 0, h, &M[26 * h * h], n);           //C12 = M[26]
		copybackM(c, h, 0, &M[27 * h * h], n);           //C21 = M[27]
		copybackM(c, h, h, &M[28 * h * h], n);           //C22 = M[28]
		delete M;
	}
}

void Writer(int** a, int n) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			cin >> a[i][j];
}

void AutoWriter(int** a, int n) {
	for (int i = 0; i < n; i++) {
		int* G = new int[ch];
		for (int j = 0; j < ch; j++)
			G[j] = 0;
		for (int j = 0; j < n; j++)
			G[j] = rand() % 10;
		a[i] = G;
	}
}

void AutoZero(int** a, int n) {
	for (int i = 0; i < n; i++) {
		int* G = new int[n];
		for (int j = 0; j < n; j++)
			G[j] = 0;
		a[i] = G;
	}
}

void Reader(int** a, int n) {
	for (int i = 0; i < n; i++) {
		int* G = a[i];
		cout << endl;
		for (int j = 0; j < n; j++)
			cout << G[j] << " ";
	}
	cout << endl;
}

int main() {
	setlocale(LC_ALL, "Russian");
	
	int** a;
	int** b;
	int** c;
	int** c1;
	int n, control;

	cout << "Введите размерность массива:" << endl;
	cin >> n;

	while (ch < n)
		ch = ch * 2;

	a = new int* [ch];
	for (int i = 0; i < ch; i++)
		a[i] = new int[ch];
	AutoZero(a, ch);

	b = new int* [ch];
	for (int i = 0; i < ch; i++)
		b[i] = new int[ch];
	AutoZero(b, ch);

	c = new int* [ch];
	for (int i = 0; i < ch; i++)
		c[i] = new int[ch];
	AutoZero(c, ch);

	c1 = new int* [ch];
	for (int i = 0; i < ch; i++)
		c1[i] = new int[ch];
	AutoZero(c1, ch);

	cout << endl << "Каким способом хотите заполнить элементы массивов:" << endl << "(1) В ручную" << endl << "(2) Автоматически" << endl;
	cin >> control;
	if (control == 1) {
		cout << "Заполните массив А:" << endl;
		Writer(a, n);
		cout << "Заполните массив В:" << endl;
		Writer(b, n);
		cout << "Массив А:" << endl;
		Reader(a, n);
		cout << "Массив B:" << endl;
		Reader(b, n);
	} else {
		AutoWriter(a, n);
		AutoWriter(b, n);
		cout << "Массив А:" << endl;
		Reader(a, ch);
		cout << "Массив B:" << endl;
		Reader(b, ch);
	}
	unsigned int start_time1 = clock();
	mul_normal(c1, a, b, n);
	unsigned int end_time1 = clock();
	unsigned int search_time1 = end_time1 - start_time1;

	cout << "Произведение A*B = C(Тривиальным методом):" << endl;
	Reader(c1, n);
	unsigned int start_time2 = clock();
	mul(c, a, b, ch);
	unsigned int end_time2 = clock();
	unsigned int search_time2 = end_time2 - start_time2;

	cout << "Произведение A*B = C(методом Штрассена):" << endl;
	Reader(c, n);
	cout << "Время работы тривиального метода = " << search_time1 << endl;
	cout << "Время работы метода Штрассена = " << search_time2 << endl;

	for (int i = 0; i < ch; i++)
		delete[] a[i];
	delete[] a;
	for (int i = 0; i < ch; i++)
		delete[] b[i];
	delete[] b;
	for (int i = 0; i < ch; i++)
		delete[] c[i];
	delete[] c;
	for (int i = 0; i < ch; i++)
		delete[] c1[i];
	delete[] c1;
	return 0;
}