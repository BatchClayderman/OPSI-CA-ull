#include <iostream>
#include <vector>
#include <algorithm>
#ifndef _OPSICA_H
#define _OPSICA_H
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EOF
#define EOF (-1)
#endif
#ifndef NULL
#define NULL 0
#endif
#ifndef MAX_BUFFER
#define MAX_BUFFER 5000
#endif
#ifndef MAX_PATH
#ifdef _MAX_PATH
#define MAX_PATH _MAX_PATH
#else
#define MAX_PATH 260
#endif
#endif

#define kBit 128
#define N 26
#define n 11
#define beta n
#define gamma 3
#define sigma 40
#define lambda 128
#define sizeW 128
#define TimeToTest 50
#endif//_OPSICA_H
using namespace std;
typedef unsigned long long int Element;
typedef const void* CPVOID;
vector<int> hashpi{}, archashpi{};
size_t baseNum = kBit / (sizeof(Element) << 3);
clock_t sub_start_time = clock(), sub_end_time = clock();
double timerR = 0, timerS = 0, timerC = 0;


/* 子函数 */
void init_hashpi()
{
	for (int i = 0; i < beta; ++i)
	{
		hashpi.push_back(i);
		archashpi.push_back(NULL);// initial
	}
	random_shuffle(hashpi.begin(), hashpi.end());
	for (int i = 0; i < beta; ++i)
		archashpi[hashpi[i]] = i;
	return;
}

int pi(int index)
{
	return hashpi[index % beta];
}

int arcpi(int value)
{
	return archashpi[value % beta];
}

Element getRandom()//获取随机数
{
	Element random = rand();
	random <<= 32;
	random += rand();
	return random;
}

void randBool(bool** arr, size_t a, size_t b)
{
	time_t t;
	srand((unsigned int)time(&t));
	for (size_t i = 0; i < a; ++i)
		for (size_t j = 0; j < b; ++j)
			*((bool*)arr + i * a + j) = (bool)(rand() % 2);
	return;
}

Element F(Element a, Element b)
{
	return a ^ b;// example
}

void getInput(Element array[], int size)//获得输入
{
	char buffer[MAX_BUFFER] = { 0 }, cTmp[MAX_PATH] = { 0 };
	rewind(stdin);
	fflush(stdin);
	fgets(buffer, MAX_BUFFER, stdin);
	int cIndex = 0, eIndex = 0;
	for (int i = 0; i < MAX_BUFFER; ++i)
		if (buffer[i] >= '0' && buffer[i] <= '9')
			cTmp[cIndex++] = buffer[i];
		else if (cIndex)
		{
			char* endPtr;
			array[eIndex++] = strtoull(cTmp, &endPtr, 0);
			if (eIndex >= size)
				return;
			cIndex = 0;// Rewind cIndex
			memset(cTmp, 0, strlen(cTmp));// Rewind cTmp
		}
	return;
}

Element r_i(Element ele, int i)//哈希函数
{
	return ele << i;
}

int compare(CPVOID a, CPVOID b)//比较函数
{
	return (int)(*(Element*)a - *(Element*)b);
}

int BinarySearch(Element array_lists[], int nBegin, int nEnd, Element target, unsigned int& compareCount)
{
	if (nBegin > nEnd)
		return EOF;//未能找到目标
	int nMid = (nBegin + nEnd) >> 1;//使用位运算加速
	++compareCount;
	if (array_lists[nMid] == target)//找到目标
		return nMid;
	else if (array_lists[nMid] > target)//分而治之
		return BinarySearch(array_lists, nBegin, nMid - 1, target, compareCount);
	else
		return BinarySearch(array_lists, nMid + 1, nEnd, target, compareCount);
}

Element factorial(Element num)
{
	Element result = 1;
	for (Element i = 1; i <= num; ++i)
		result *= i;
	return result;
}

Element combination(Element a, Element b)
{
	return factorial(a) / (factorial(b) * factorial(a - b));
}

double Cp1p(Element a, Element b, double p)
{
	return (double)combination(a, b) * pow(p, b) * pow(1 - p, a - b);
}

#ifdef _DEBUG
void printArray(Element* arr, size_t length, string name)
{
	if (length <= 0 || NULL == arr)
	{
		cout << name << " = {}" << endl;
		return;
	}
	cout << name << " = { " << *arr;
	for (int i = 1; i < length; ++i)
		cout << ", " << arr[i];
	cout << " }" << endl;
	return;
}
#else
void printArray()
{
	return;
}
#endif


/* 类 */
class Receiver
{
private:
	Element X[n] = { NULL };
	Element X_c[beta] = { NULL };
	Element Z[beta] = { NULL };
	Element Z_pi[beta] = { NULL };
	Element V[beta] = { NULL };
	Element W[beta] = { NULL };
	vector<Element> intersection{};

public:
	void input_X()
	{
		getInput(this->X, n);
		return;
	}
	void auto_input_X()
	{
		vector<Element> v;
		while (v.size() < n)
		{
			Element tmp = getRandom();
			if (find(v.begin(), v.end(), tmp) == v.end())
				v.push_back(tmp);
		}
		for (int i = 0; i < n; ++i)
			this->X[i] = v[i];
		return;
	}
	void hash_X_to_X_c()
	{
		for (int i = 0; i < n; ++i)
		{
			int index = r_i(this->X[i], 1) % beta;
			if (0 != this->X_c[index])// already exist
			{
				int new_index = r_i(this->X_c[index], 2) % beta;
				if (0 != this->X_c[new_index])// still already exist
					;// abundant
				else
					this->X_c[new_index] = this->X_c[index];
			}
			else
				this->X_c[index] = X[i];
		}
		return;
	}
	void splitXc()
	{
		for (int i = 0; i < n; ++i)
		{
			this->Z[i] = getRandom();
			this->Z_pi[i] = this->X_c[i] ^ this->Z[i];
		}
		return;
	}
	Element* send_Z()
	{
		return this->Z;
	}
	Element* send_Z_pi()
	{
		return this->Z_pi;
	}
	void receive_V(Element* V)
	{
		for (int i = 0; i < beta; ++i)
			this->V[i] = *(V + i);
		return;
	}
	void receive_W(Element* W)
	{
		for (int i = 0; i < beta; ++i)
			this->W[i] = *(W + i);
		return;
	}
	void printIntersection()
	{
		this->intersection.clear();
		unsigned int comparation = 0;
		qsort(this->W, beta, sizeof(Element), compare);
		qsort(this->V, beta, sizeof(Element), compare);
		for (int i = 0; i < beta; ++i)
			if (this->W[i] && BinarySearch(this->V, 0, beta - 1, this->W[i], comparation))
				this->intersection.push_back(this->W[i]);
		if (this->intersection.size())
		{
#ifdef _DEBUG
			cout << "| V ∩ W | = | { " << intersection[0];
			for (size_t i = 1; i < this->intersection.size(); ++i)
				cout << ", " << this->intersection[i];
			cout << " } | = " << this->intersection.size() << endl;
#else
			cout << "| V ∩ W | = " << this->intersection.size() << endl;
#endif
		}
		else
			cout << "| V ∩ W | = 0" << endl;
		return;
	}
	size_t printSize(bool isPow)
	{
		cout << "Timeof(R) = " << timerR * baseNum << " / " << TimeToTest << " = " << timerR * baseNum / TimeToTest << " ms" << endl;
		cout << "sizeof(Receiver) = " << sizeof(Receiver) * baseNum << " KB" << endl;
		cout << "sizeof(R) = " << sizeof(this) * baseNum << " MB" << endl;
		cout << "\tsizeof(R.X) = " << (isPow ? (sizeof(Element) * baseNum) << n : sizeof(this->X) * baseNum) << " B" << endl;
		cout << "\tsizeof(R.X_c) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->X_c) * baseNum) << " B" << endl;
		cout << "\tsizeof(R.Z) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->Z) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(R.Z_pi) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->Z_pi) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(R.V) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->V) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(R.W) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->W) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(R.intersection) = " << (isPow ? (sizeof(Element) * baseNum) << this->intersection.size() : sizeof(this->intersection) * baseNum) << " B" << endl;
		cout << "\tsizeof(R.*) = " << (isPow ? sizeof(Element) * baseNum * ((1 << beta) + (1 << beta) + (1 << beta) + (1 << beta)) : (sizeof(this->Z) + sizeof(this->Z_pi) + sizeof(this->V) + sizeof(this->W)) * baseNum) << " B (*)" << endl;
		return isPow ? sizeof(Element) * baseNum * ((1 << beta) + (1 << beta) + (1 << beta) + (1 << beta)) : (sizeof(this->Z) + sizeof(this->Z_pi) + sizeof(this->V) + sizeof(this->W)) * baseNum;
	}
};
Receiver R;

class Sender
{
private:
	Element Y[N] = { NULL };
	Element Z_pi[beta] = { NULL };
	Element I[N] = { NULL };
	bool A[beta][sizeW] = { { false } };
	Element k = NULL;
	Element V[beta] = { NULL };
	Element omega = NULL;
	Element l1 = lambda << 1, l2 = sigma + (Element)log(beta * beta);

public:
	void input_Y()
	{
		getInput(this->Y, N);
	}
	void auto_input_Y()
	{
		vector<Element> v;
		while (v.size() < N)
		{
			Element tmp = getRandom();
			if (find(v.begin(), v.end(), tmp) == v.end())
				v.push_back(tmp);
		}
		for (int i = 0; i < N; ++i)
			this->Y[i] = v[i];
		return;
	}
	void receive_Z_pi(Element* Z)
	{
		for (int i = 0; i < beta; ++i)
			this->Z_pi[i] = *(Z + i);
		return;
	}
	void compute_Q_I()
	{
		for (int i = 0; i < N; ++i)
		{
			Element q_j = 0, I_i = 0;
			for (int j = 0; j < gamma; ++j)
			{
				q_j = arcpi(r_i(this->Y[i], j) % beta);
				this->I[i] = this->Y[i] ^ this->Z_pi[q_j];
			}
		}
		return;
	}
	void* sendNULL()
	{
		randBool((bool**)this->A, beta, sizeW);
		return NULL;
	}
	void obtain_k(bool** B, Element* W)
	{
		this->k = getRandom();
		for (int i = 0; i < beta; ++i)
		{
			size_t a = (size_t)(this->k % (beta * sizeW) / beta);
			size_t b = this->k % (beta * sizeW) % beta;
			W[i] = (*((bool*)this->A + a * n + b) ^ *((bool*)B + a * n + b)) + (this->k >> (beta - i));
		}
		return;
	};
	void get_V()
	{
		for (int i = 0; i < beta; ++i)
			this->V[i] = F(this->k, this->I[i % N]);
		return;
	}
	Element* send_V()
	{
		return this->V;
	}
	void Cp1p_omega()
	{
		this->omega = NULL;
		double p = 0.5, total = 0;
		while (total <= p)
		{
			total += Cp1p(n, this->omega, p);
			++this->omega;
		}
		return;
	}
	size_t printSize(bool isPow)
	{
		cout << "Timeof(S) = " << timerS * baseNum << " / " << TimeToTest << " = " << timerS * baseNum / TimeToTest << " ms" << endl;
		cout << "sizeof(Sender) = " << sizeof(Sender) << " KB" << endl;
		cout << "sizeof(S) = " << sizeof(this) * baseNum << " MB" << endl;
		cout << "\tsizeof(S.Y) = " << (isPow ? (sizeof(Element) * baseNum) << N : sizeof(this->Y) * baseNum) << " B" << endl;
		cout << "\tsizeof(S.Z_pi) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->Z_pi) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(S.I) = " << (isPow ? (sizeof(Element) * baseNum) << N : sizeof(this->I) * baseNum) << " B" << endl;
		cout << "\tsizeof(S.A) = " << (isPow ? (sizeof(Element) * baseNum) * (beta * sizeW) : sizeof(this->A) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(S.k) = " << sizeof(this->k) * baseNum << " B (*)" << endl;
		cout << "\tsizeof(S.V) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->V) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(S.*) = " << (isPow ? (sizeof(Element) * baseNum) * ((1 << beta) + (beta * sizeW) + 1 + (1 << beta)) : (sizeof(this->Z_pi) + sizeof(this->A) + sizeof(this->k) + sizeof(this->V)) * baseNum) << " B (*)" << endl;
		return isPow ? (sizeof(Element) * baseNum) * ((1 << beta) + (beta * sizeW) + 1 + (1 << beta)) : (sizeof(this->Z_pi) + sizeof(this->A) + sizeof(this->k) + sizeof(this->V)) * baseNum;
	}
};
Sender S;

class Cloud
{
private:
	Element Z[beta] = { NULL };
	bool B[beta][sizeW] = { { false } };
	Element W[beta] = { NULL };
	Element l1 = lambda << 1, l2 = sigma + (Element)log(beta * beta);

public:
	void receive_Z(Element* Z)
	{
		for (int i = 0; i < beta; ++i)
			this->Z[i] = *(Z + i);
		return;
	}
	void receiveNULL(void* pNull)
	{
		pNull;
		return;
	}
	bool** act_Z()
	{
		randBool((bool**)this->B, beta, sizeW);
		return (bool**)this->B;
	}
	Element* obtain_W()
	{
		return this->W;
	}
	Element* send_W()
	{
		return this->W;
	}
	size_t printSize(bool isPow)
	{
		cout << "Timeof(C) = " << timerC * baseNum << " / " << TimeToTest << " = " << timerC * baseNum / TimeToTest << " ms" << endl;
		cout << "sizeof(Cloud) = " << sizeof(Cloud) << " KB" << endl;
		cout << "sizeof(C) = " << sizeof(this) * baseNum << " MB" << endl;
		cout << "\tsizeof(C.Z) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->Z) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(C.B) = " << (isPow ? (sizeof(Element) * baseNum) * (beta * sizeW) : sizeof(this->B) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(C.W) = " << (isPow ? (sizeof(Element) * baseNum) << beta : sizeof(this->W) * baseNum) << " B (*)" << endl;
		cout << "\tsizeof(C.*) = " << (isPow ? (sizeof(Element) * baseNum) * ((1 << beta) + (beta * sizeW) + (1 << beta)) : (sizeof(this->Z) + sizeof(this->B) + sizeof(this->W)) * baseNum) << " B (*)" << endl;
		return isPow ? (sizeof(Element) * baseNum) * ((1 << beta) + (beta * sizeW) + (1 << beta)) : (sizeof(this->Z) + sizeof(this->B) + sizeof(this->W)) * baseNum;
	}
};
Cloud C;


/* 主要函数 */
void initial(bool isAuto)
{
	init_hashpi();// setup pi
	if (isAuto)
	{
		S.auto_input_Y();
		R.auto_input_X();
	}
	else
	{
		cout << "Please input array Y with size " << N << ": " << endl;
		S.input_Y();// Sender S has input Y
		cout << endl;
		cout << "Please input array X with size " << n << ": " << endl;
		R.input_X();// Sender R has input X
		cout << endl;
	}
	return;
}

void protocol()
{
	/* Step 1 */
	sub_start_time = clock();
	R.hash_X_to_X_c();
	sub_end_time = clock();
	timerR += (double)sub_end_time - sub_start_time;
	
	/* Step 2 */
	sub_start_time = clock();
	R.splitXc();
	sub_end_time = clock();
	timerR += (double)sub_end_time - sub_start_time;
	
	/* Step 3 */
	sub_start_time = clock();
	C.receive_Z(R.send_Z());
	sub_end_time = clock();
	timerR += (double)sub_end_time - sub_start_time;
	timerC += (double)sub_end_time - sub_start_time;// Z is sent to C
	sub_start_time = clock();
	S.receive_Z_pi(R.send_Z_pi());
	sub_end_time = clock();
	timerR += (double)sub_end_time - sub_start_time;
	timerS += (double)sub_end_time - sub_start_time;// Z_pi is sent to C
	
	/* Step 4 */
	sub_start_time = clock();
	S.compute_Q_I();
	sub_end_time = clock();
	timerS += (double)sub_end_time - sub_start_time;
	
	/* Step 5 */
	sub_start_time = clock();
	C.receiveNULL(S.sendNULL());
	S.obtain_k(C.act_Z(), C.obtain_W());
	sub_end_time = clock();
	timerS += (double)sub_end_time - sub_start_time;
	timerC += (double)sub_end_time - sub_start_time;
	
	/* Step 6 */
	sub_start_time = clock();
	S.get_V();
	S.Cp1p_omega();
	sub_end_time = clock();
	timerS += (double)sub_end_time - sub_start_time;
	sub_start_time = clock();
	R.receive_V(S.send_V());
	sub_end_time = clock();
	timerS += (double)sub_end_time - sub_start_time;
	timerR += (double)sub_end_time - sub_start_time;

	/* Step 7 */
	sub_start_time = clock();
	R.receive_W(C.send_W());
	sub_end_time = clock();
	timerR += (double)sub_end_time - sub_start_time;
	timerC += (double)sub_end_time - sub_start_time;

	/* Step 8 */
	sub_start_time = clock();
	R.printIntersection();
	sub_end_time = clock();
	timerR += (double)sub_end_time - sub_start_time;
	return;
}



/* main 函数 */
void test()
{
	time_t t;
	srand((unsigned int)time(&t));
	initial(true);
	protocol();
	return;
}

int main()
{
	clock_t start_time = clock();
	for (int i = 0; i < TimeToTest; ++i)
	{
		cout << "/**************************************** Round: " << i + 1 << " ****************************************/" << endl;
		test();
		cout << endl << endl;
	}
	clock_t end_time = clock();
	cout << endl;
	cout << "/**************************************** OPSI-CA ****************************************/" << endl;
	cout << "kBit = " << kBit << "\t\tgamma = " << gamma << "\t\tlambda = " << lambda << endl;
	cout << "beta = [2 ** " << (log2(1.27) + beta) << "]\t\tN = 2 ** " << N << "\t\tn = 2 ** " << n << endl;
	cout << "Time: " << ((double)end_time - start_time) * baseNum << " / " << TimeToTest << " = " << ((double)end_time - start_time) * baseNum / TimeToTest << "ms" << endl;
	cout << "sizeof(*) = " << (R.printSize(true) + S.printSize(true) + C.printSize(true)) << " B (*)" << endl << endl;
	return EXIT_SUCCESS;
}