// EASYNN.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <map>
#include <string>
#include <sstream>
using namespace std;

class doublevector {
	int _size;
	int max_size;
	int * count;
public:
	double * data;
	doublevector() : data(NULL), _size(0), max_size(0), count(NULL) {}
	doublevector(int n) :_size(n), max_size(n) {
		data = new double[n];
		memset(data, 0, n*sizeof(double));
		count = new int;
		*count = 1;
	}
	doublevector(const doublevector & other) : _size(other.size()), max_size(other.max_size), count(other.count) {
		data = other.data;
		++(*count);
	}
	doublevector operator = (const doublevector & other) {

		if (data != NULL && (--(*count)) == 0) {
			delete[] data;
			delete count;
		}

		_size = other.size();
		max_size = other.max_size;
		data = other.data;
		count = other.count;
		++(*count);
		return *this;
	}
	doublevector copy() {
		doublevector ans;
		ans.data = new double[max_size];
		memcpy(ans.data, data, _size * sizeof(double));
		ans.count = new int;
		*ans.count = 1;
		return ans;
	}
	doublevector(int n, double _d) :_size(n), max_size(n) {
		data = new double[n];
		for (int i = 0; i < n; ++i)
			data[i] = _d;
		count = new int;
		*count = 1;
	}
	double & operator [] (int n) {
		return data[n];
	}
	int size() const {
		return _size;
	}
	void clear() {
		_size = 0;
	}
	bool empty() {
		return _size == 0;
	}
	void push_back(double _d) {
		if (_size < max_size)
			data[_size++] = _d;
		else {
			max_size = 2 * max_size + 2;
			double * newdata = new double[max_size];
			memcpy(newdata, data, _size * sizeof(double));
			if (data != NULL && --(*count) == 0)
				delete[] data;
			else
				count = new int;
			*count = 1;
			data = newdata;
			data[_size++] = _d;
		}
	}
	void resize(int n) {
		if (max_size >= n) {
			_size = n;
		}
		else {
			max_size = 2 * n + 2;
			double * newdata = new double[max_size];
			memcpy(newdata, data, _size * sizeof(double));
			if (--(*count) == 0)
				delete[] data;
			else
				count = new int;
			*count = 1;
			data = newdata;
			_size = n;
		}
	}
	~doublevector() {
		if (data == NULL)
			return;
		if ((--(*count)) == 0) {
			delete[] data;
			delete count;
		}
		return;
	}
};

class cutility {
public:
	static double sum(doublevector & _v) {
		double s = 0;
		for (unsigned int i = 0; i < _v.size(); ++i) {
			s += _v[i];
		}
		return s;
	}
	static doublevector add(doublevector & _v1, doublevector & _v2) {
		if (_v1.size() == 0)
			return _v2;
		if (_v1.size() != _v2.size()) {
			if (_v1.size() == 1)
				return add(doublevector(_v2.size(), _v1[0]), _v2);
			else if (_v2.size() == 1)
				return add(_v1, doublevector(_v1.size(), _v2[0]));
			throw(string("utility.add error"));
		}
		doublevector ans(_v1.size());
		for (unsigned int i = 0; i < _v1.size(); ++i) {
			ans[i] = _v1[i] + _v2[i];
		}
		return ans;
	}
	static doublevector subtract(doublevector & _v1, doublevector & _v2) {
		if (_v1.size() != _v2.size()) {
			if (_v1.size() == 1)
				return subtract(doublevector(_v2.size(), _v1[0]), _v2);
			else if (_v2.size() == 1)
				return subtract(_v1, doublevector(_v1.size(), _v2[0]));
			throw(string("utility.add error"));
		}
		doublevector ans(_v1.size());
		for (unsigned int i = 0; i < _v1.size(); ++i) {
			ans[i] = _v1[i] - _v2[i];
		}
		return ans;
	}
	static doublevector subtract(doublevector & _v1) {
		return cutility::subtract(doublevector(1, 0), _v1);
	}
	static doublevector multiply(doublevector & _v1, doublevector & _v2) {
		if (_v1.size() == 0)
			return _v2;
		if (_v1.size() != _v2.size()) {
			if (_v1.size() == 1)
				return multiply(doublevector(_v2.size(), _v1[0]), _v2);
			else if (_v2.size() == 1)
				return multiply(_v1, doublevector(_v1.size(), _v2[0]));
			throw(string("utility.add error"));
		}
		doublevector ans(_v1.size());
		for (unsigned int i = 0; i < _v1.size(); ++i) {
			ans[i] = _v1[i] * _v2[i];
		}
		return ans;
	}
	static doublevector abs(doublevector & _v1) {
		doublevector ans(_v1.size());
		for (unsigned int i = 0; i < _v1.size(); ++i) {
			ans[i] = fabs(_v1[i]);
		}
		return ans;
	}
	static double gaussrand()
	{
		static double V1, V2, S;
		static int phase = 0;
		double X;

		if (phase == 0) {
			do {
				double U1 = (double)rand() / RAND_MAX;
				double U2 = (double)rand() / RAND_MAX;

				V1 = 2 * U1 - 1;
				V2 = 2 * U2 - 1;
				S = V1 * V1 + V2 * V2;
			} while (S >= 1 || S == 0);

			X = V1 * sqrt(-2 * log(S) / S);
		}
		else
			X = V2 * sqrt(-2 * log(S) / S);

		phase = 1 - phase;

		return X;
	}
};

class node;
class NN;

class connector {
public:
	node * pre, *succ;
	doublevector forward_value;
	doublevector backward_value;
	vector<doublevector > packed_forward_value;
	vector<doublevector > packed_backward_value;
	void clear() {
		forward_value.clear();
		backward_value.clear();
	}
};

class node {
public:
	static const int START = 0;
	static const int ADD = 1;
	static const int MULTI = 2;
	static const int SUB = 3;
	static const int CONST = 4;
	static const int END = 5;
	static const int VAR = 6;
	static const int ABS = 7;
	static const int INPUT = 8;
	static const int PACK = 9;
	static const int UNPACK = 10;
	static const int COMPLEX = 11;
	static const int ACTIVE = 12;
	static int currentid;
	string name;
	unsigned int type;


	string left;
	string right;
	string middle;
	node(int _type) : type(_type) {
		switch (type) {
		case START:
			break;
		case ADD:
			middle = "+";
			break;
		case SUB:
			middle = "-";
			break;
		case MULTI:
			middle = "*";
			break;
		case ABS:
			left = "|";
			right = "|";
			break;
		}
		id = currentid++;
	}
	node(int _type, node * p);
	node(int _type, node * p, node * q);
	node(int _type, vector<node *> _v);

	unsigned int id;

	vector<connector *> pre_nodes;
	vector<connector *> succ_nodes;
	void clear() {
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			pre_nodes[i]->clear();
		}
		for (unsigned int i = 0; i < succ_nodes.size(); ++i) {
			succ_nodes[i]->clear();
		}
	}
	virtual bool forward() = 0;
	virtual bool backward() = 0;

	doublevector accumulate_backward() {
		doublevector ans;
		for (unsigned int i = 0; i < succ_nodes.size(); ++i) {
			connector *p = succ_nodes[i];
			ans = cutility::add(ans, p->backward_value);
		}
		return ans;
	}
	vector<doublevector > accumulate_packed_backward() {
		vector<doublevector > ans;
		for (int j = 0; j < succ_nodes[0]->packed_backward_value.size(); ++j) {
			doublevector temp;
			for (unsigned int i = 0; i < succ_nodes.size(); ++i) {
				connector *p = succ_nodes[i];
				temp = cutility::add(temp, p->packed_backward_value[j]);
			}
			ans.push_back(temp);
		}
		return ans;
	}
	void allocate_forward(doublevector & _forward_value) {
		for (unsigned int i = 0; i < succ_nodes.size(); ++i) {
			connector *p = succ_nodes[i];
			p->forward_value = _forward_value;
		}
	}
	void allocate_packed_forward(vector<doublevector > _packed_forward_value) {
		for (unsigned int i = 0; i < succ_nodes.size(); ++i) {
			connector *p = succ_nodes[i];
			p->packed_forward_value = _packed_forward_value;
		}
	}
	virtual void SGD(double alpha) {
		return;
	}
	virtual void momentumSGD(double lambda, double alpha) {
		return;
	}
	~node() {
		for (int i = 0; i < succ_nodes.size(); ++i) {
			connector *p = succ_nodes[i];
			while (true) {
				vector<connector *>::iterator iter = find(p->succ->pre_nodes.begin(), p->succ->pre_nodes.end(), p);
				if (iter == p->succ->pre_nodes.end())
					break;
				p->succ->pre_nodes.erase(iter);
			}
			delete p;
		}
	}
	virtual void init(void * data = NULL) {
		return;
	}
	virtual void absorb(node * other) {
		return;
	}

	virtual node * copy() = 0;
};

int node::currentid = 0;



class startnode : public node {
public:
	startnode() :node(START) {}
	map<int, doublevector > data;
	bool forward() {
		for (unsigned int i = 0; i < succ_nodes.size(); ++i) {
			if (data.find(succ_nodes[i]->succ->id) != data.end()) {
				succ_nodes[i]->forward_value = data[succ_nodes[i]->succ->id];
			}
		}
		return true;
	}
	bool backward() {
		return true;
	}
	virtual startnode * copy() {
		return new startnode();
	}
};

class addnode : public node {
public:
	addnode() :node(ADD) {}
	addnode(node * p, node * q) : node(ADD, p, q) {}
	addnode(node * p) : node(ADD, p) {}
	addnode(vector<node *> v) : node(ADD, v) {}
	bool forward() {
		doublevector forward_value;
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			forward_value = cutility::add(forward_value, p->forward_value);
		}
		allocate_forward(forward_value);
		return !forward_value.empty();
	}
	bool backward() {
		double backward_value = 0;
		doublevector info = accumulate_backward();
		bool has_backward = !info.empty();
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			p->backward_value = info;
		}
		return has_backward;
	}
	virtual addnode * copy() {
		return new addnode();
	}
};

class multiplynode : public node {
public:
	multiplynode() :node(MULTI) {}
	multiplynode(node * p, node * q) : node(MULTI, p, q) {}
	multiplynode(node * p) : node(MULTI, p) {}
	multiplynode(vector<node *> v) : node(MULTI, v) {}
	bool forward() {
		doublevector forward_value;
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			forward_value = cutility::multiply(forward_value, p->forward_value);
		}
		allocate_forward(forward_value);
		return forward_value.size() > 0;
	}
	bool backward() {
		double backward_value = 0;
		doublevector info = accumulate_backward();
		bool has_backward = info.size() > 0;
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			doublevector temp_value;
			for (int j = 0; j < pre_nodes.size(); ++j) {
				if (j != i) {
					temp_value = cutility::multiply(temp_value, pre_nodes[j]->forward_value);
				}
			}
			p->backward_value = cutility::multiply(info, temp_value);
		}
		return has_backward;
	}
	virtual multiplynode * copy() {
		return new multiplynode();
	}
};

class constantnode : public node {
public:
	constantnode() :node(CONST) {}
	int value;
	bool forward() {
		allocate_forward(doublevector(1, value));
		return true;
	}
	bool backward() {
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			p->backward_value = doublevector();
		}
		return true;
	}
	virtual constantnode * copy() {
		constantnode * ans = new constantnode();
		ans->value = value;
		return ans;
	}
};

class subtractnode : public node {
public:
	subtractnode() :node(SUB) {}
	subtractnode(node * p, node * q) : node(SUB, p, q) {}
	subtractnode(node * p) : node(SUB, p) {}
	bool forward() {
		doublevector forward_value;
		if (pre_nodes.size() == 1) {
			forward_value = cutility::subtract(pre_nodes[0]->forward_value);
		}
		else if (pre_nodes.size() == 2) {
			forward_value = cutility::subtract(pre_nodes[0]->forward_value, pre_nodes[1]->forward_value);
		}
		else {
			throw(string("subtractnode: 2 inputs needed"));
		}
		allocate_forward(forward_value);
		return forward_value.size() > 0;
	}
	bool backward() {
		doublevector info = accumulate_backward();
		if (pre_nodes.size() == 1) {
			pre_nodes[0]->backward_value = cutility::subtract(info);
		}
		else if (pre_nodes.size() == 2) {
			pre_nodes[0]->backward_value = info;
			pre_nodes[1]->backward_value = cutility::subtract(info);
		}
		else {
			throw(string("subtractnode: 2 inputs needed"));
		}
		return info.size() > 0;
	}
	virtual subtractnode * copy() {
		return new subtractnode();
	}
};

class abstractnode : public node {
public:
	abstractnode() :node(ABS) {}
	abstractnode(node * p) : node(ABS, p) {}
	bool forward() {
		doublevector forward_value;
		if (pre_nodes.size() == 1) {
			forward_value = cutility::abs(pre_nodes[0]->forward_value);
		}
		else {
			throw(string("absnode: 1 input needed"));
		}
		allocate_forward(forward_value);
		return forward_value.size() > 0;
	}
	bool backward() {
		doublevector info = accumulate_backward();
		if (pre_nodes.size() == 1) {
			if (info.size() == 1) {
				info = doublevector(pre_nodes[0]->forward_value.size(), info[0]);
			}
			if (pre_nodes[0]->forward_value.size() != info.size()) {
				pre_nodes[0]->backward_value.clear();
				return false;
			}
			pre_nodes[0]->backward_value.resize(info.size());
			for (unsigned int i = 0; i < info.size(); ++i) {
				pre_nodes[0]->backward_value[i] = (pre_nodes[0]->forward_value[i] >= 0 ? info[i] : -info[i]);
			}
		}
		else {
			throw(string("absnode: 1 inputs needed"));
		}
		return true;
	}
	virtual abstractnode * copy() {
		return new abstractnode();
	}
};

class endnode : public node {
public:
	doublevector error;
	double averageerror;
	endnode() :node(END) {}
	endnode(node * p) : node(END, p) {
		type = END;
	}
	bool forward() {
		doublevector forward_value;
		if (pre_nodes.size() == 1) {
			error = pre_nodes[0]->forward_value;
		}
		else {
			throw(string("absnode: 1 input needed"));
		}
		averageerror = cutility::sum(error) / error.size();
		return forward_value.size() > 0;
	}
	bool backward() {
		if (pre_nodes.size() == 1) {
			pre_nodes[0]->backward_value = doublevector(1, 1);
		}
		else {
			throw(string("absnode: 1 inputs needed"));
		}
		return !error.empty();
	}
	virtual endnode * copy() {
		return new endnode();
	}
};

class inputnode : public addnode {
public:
	inputnode() {
		type = INPUT;
	}
	inputnode(node * p) : addnode(p) {
		type = INPUT;
	}
	virtual inputnode * copy() {
		return new inputnode();
	}
};

class variablenode : public node {
public:
	variablenode() :node(VAR) {
		velocity = 0;
		diff = 0;
		value = 0;
		type = VAR;
	}
	variablenode(node * p) : node(VAR, p) {
		velocity = 0;
		diff = 0;
		value = 0;
		type = VAR;
	}
	double value;
	double diff;
	double velocity;
	bool forward() {
		allocate_forward(doublevector(1, value));
		return true;
	}
	bool backward() {
		double backward_value = 0;
		doublevector info = accumulate_backward();
		if (info.empty())
			return false;
		diff = cutility::sum(info) / info.size();
		return true;
	}
	virtual void SGD(double alpha) {
		velocity = alpha * diff;
		value = value - alpha * diff;
	}
	virtual void momentumSGD(double lambda, double alpha) {
		velocity = lambda * velocity - alpha * diff;
		value = value + velocity;
	}
	virtual void init(void * data = NULL) {
		if (data == NULL)
			value = cutility::gaussrand();
		else {
			value = (cutility::gaussrand() * ((double *)data)[1]) + ((double *)data)[0];
		}
	}
	virtual void absorb(node * other) {
		if (other->type != node::VAR)
			return;
		for (int i = 0; i < other->succ_nodes.size(); ++i) {
			for (int j = 0; j < other->succ_nodes[i]->succ->pre_nodes.size(); ++j) {
				if (other->succ_nodes[i]->succ->pre_nodes[j]->pre == other) {
					other->succ_nodes[i]->succ->pre_nodes[j]->pre = this;
				}
			}
			succ_nodes.push_back(other->succ_nodes[i]);
		}
		delete other;
	}
	virtual variablenode * copy() {
		return new variablenode();
	}
};

class packnode : public node {
public:
	packnode() : node(PACK) {}
	bool forward() {
		vector<doublevector > packed_forward_value;
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			packed_forward_value.push_back(p->forward_value);
		}
		allocate_packed_forward(packed_forward_value);
		return !packed_forward_value.empty();
	}
	bool backward() {
		vector<doublevector > info = accumulate_packed_backward();
		if (info.size() != pre_nodes.size())
			throw(string("packnode error"));
		for (unsigned int i = 0; i < pre_nodes.size(); ++i) {
			connector * p = pre_nodes[i];
			p->backward_value = info[i];
		}
		return true;
	}
	virtual packnode * copy() {
		return new packnode();
	}
};

class unpacknode : public node {
public:
	unpacknode() : node(PACK) {}
	bool forward() {
		if (pre_nodes.size() != 1)
			throw(string("unpacknode: more than one input"));
		if (pre_nodes[0]->packed_forward_value.size() != succ_nodes.size())
			throw(string("unpacknode: number error"));
		for (int i = 0; i < succ_nodes.size(); ++i) {
			succ_nodes[i]->forward_value = pre_nodes[0]->packed_forward_value[i];
		}
		return true;
	}
	bool backward() {
		if (pre_nodes.size() != 1)
			throw(string("unpacknode: more than one input"));
		if (pre_nodes[0]->packed_forward_value.size() != succ_nodes.size())
			throw(string("unpacknode: number error"));
		pre_nodes[0]->packed_backward_value.clear();
		for (int i = 0; i < succ_nodes.size(); ++i) {
			pre_nodes[0]->packed_backward_value.push_back(succ_nodes[i]->backward_value);
		}
		return true;
	}
	virtual unpacknode * copy() {
		return new unpacknode();
	}
};

class complexnode : public node {
public:
	complexnode() :node(COMPLEX) {}
	complexnode(node * s);
	vector<node *> v;
	void build(node * s);
	bool forward() {
		if (pre_nodes.size() != 1) {
			throw(string("complex error"));
		}
		v[0]->pre_nodes = pre_nodes;
		v[v.size() - 1]->succ_nodes = succ_nodes;
		for (node * p : v) {
			p->forward();
		}
	}
	void clear() {
		node::clear();
		for (node * p : v)
			p->clear();
		return;
	}
	bool backward() {
		if (pre_nodes.size() != 1) {
			throw(string("complex error"));
		}
		v[0]->pre_nodes = pre_nodes;
		v[v.size() - 1]->succ_nodes = succ_nodes;
		for (int i = v.size() - 1; i >= 0; --i)
			v[i]->backward();
	}
	virtual void SGD(double alpha) {
		for (node * p : v)
			p->SGD(alpha);
		return;
	}
	virtual void momentumSGD(double lambda, double alpha) {
		for (node * p : v)
			p->momentumSGD(lambda, alpha);
		return;
	}
	virtual void init() {
		for (node * p : v) {
			p->init();
		}
	}
	virtual void absorb(node * other) {
		if (other->type != COMPLEX)
			return;
		complexnode * cn = (complexnode *)other;
		for (int i = 0; i < v.size(); ++i) {
			v[i]->absorb(cn->v[i]);
		}
	}
	virtual complexnode * copy();
};

class activefunction : public node {
public:
	double(*f) (double);
	double(*fprime) (double);
	activefunction(double(*_f)(double), double(*_fprime)(double)) :node(ACTIVE), f(_f), fprime(_fprime) {
		name = "active function";
	}
	virtual bool forward() {

		doublevector x = pre_nodes[0]->forward_value;
		for (int i = 0; i < x.size(); ++i) {
			x[i] = f(x[i]);
		}
		allocate_forward(x);
		return true;
	}

	virtual bool backward() {
		doublevector dy = accumulate_backward();
		for (int i = 0; i < dy.size(); ++i) {
			dy[i] = dy[i] * fprime(pre_nodes[0]->forward_value[i]);
		}
		pre_nodes[0]->backward_value = dy;
		return true;
	}
	virtual activefunction * copy() {
		activefunction * ans = new activefunction(f, fprime);
		ans->name = name;
		return ans;
	}
};

class parallelfunction : public complexnode {
public:
	int n;
	node * template_node;
	parallelfunction(int _n, node * _t_n);
};

double FUN_RELU(double x) {
	if (x >= 0)
		return x;
	else
		return 0;
}

double FUN_RELUPRIME(double x) {
	if (x >= 0)
		return 1;
	else
		return 0;
}

double FUN_SIGMOID(double x) {
	const double E = 2.718281828459;
	double temp = pow(E, x);
	return temp / (1 + temp);
}

double FUN_SIGMOIDPRIME(double x) {
	double temp = FUN_SIGMOID(x);
	return temp*(1 - temp);
}

class RELU : public activefunction {
public:
	RELU() :activefunction(FUN_RELU, FUN_RELUPRIME) {
		name = "RELU";
	}
	virtual RELU * copy() {
		return new RELU;
	}
};

RELU * relu = new RELU;

class SIGMOID :public activefunction {
public:
	SIGMOID() : activefunction(FUN_SIGMOID, FUN_SIGMOIDPRIME) {
		name = "SIGMOID";
	}
	bool backward() {
		doublevector dy = accumulate_backward();
		for (int i = 0; i < dy.size(); ++i) {
			double temp = succ_nodes[0]->forward_value[i];
			dy[i] = dy[i] * (temp * (1 - temp));
		}
		pre_nodes[0]->backward_value = dy;
		return true;
	}
	virtual SIGMOID * copy() {
		return new SIGMOID;
	}
};

SIGMOID * sigmoid = new SIGMOID;

class LinearTransform : public complexnode {
public:
	variablenode * a, *b;
	addnode * x;
	addnode * y;
	int n, m;
	LinearTransform() {}
	LinearTransform(int n, int m);
	void test() {
		int tt = 1;
		for (node * p : v) {
			if (p->type == VAR) {
				((variablenode *)p)->value = tt++;
			}
		}
	}
	void explify() {
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				cout << a[n*j + i].value << "*";
				cout << x[i].pre_nodes[0]->forward_value[0] << "+";
			}
			cout << b[j].value << "=" << y[j].succ_nodes[0]->forward_value[0] << endl;
		}
	}
};

class NN {
public:

	static addnode * add(node * p, node * q);
	static multiplynode * multi(node * p, node * q);
	static abstractnode * abs(node * p);
	static subtractnode * sub(node * p);
	static subtractnode * sub(node * p, node * q);
	static void topologysort(node * p, vector<node *> & v);
	static void findallnodes(node * p, vector<node *> & v);
public:
	static vector<node *> topologysort(node * p);
	startnode * start;
	vector<node *> nodes;
	NN();
	static void link(node * p, node * q);
	void topologysort() {
		nodes = topologysort(start);
	}
	void run() {
		topologysort();

		int n = 10000;

		while (n--) {
			for (int i = 0; i < nodes.size(); ++i)
				nodes[i]->forward();

			for (int i = nodes.size() - 1; i >= 0; --i)
				nodes[i]->backward();

			for (int i = 0; i < nodes.size(); ++i) {
				if (nodes[i]->type == node::VAR)
					((variablenode *)nodes[i])->momentumSGD(0.1, 0.01);
			}

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->clear();
			}

		}

	}
	void test1() {
		inputnode * y = new inputnode(start);
		inputnode * x = new inputnode(start);
		inputnode * z = new inputnode(start);
		x->name = "x";
		y->name = "y";
		variablenode * a = new variablenode();
		a->init();
		a->name = "a";
		variablenode * a1 = new variablenode();
		a1->init();
		a1->name = "a1";
		doublevector xdata, ydata, zdata;
		xdata.push_back(1);
		xdata.push_back(2);
		ydata.push_back(3);
		ydata.push_back(4);
		zdata.push_back(2);
		zdata.push_back(3);
		start->data[x->id] = xdata;
		start->data[y->id] = ydata;
		start->data[z->id] = zdata;
		node * t = sub(z, add(multi(x, a), multi(y, a1)));
		a->absorb(a1); a1 = a;
		endnode * end = new endnode(multi(t, t));


		topologysort();
		int n = 3000;
		while (n--) {
			for (int i = 0; i < nodes.size(); ++i)
				nodes[i]->forward();

			for (int i = nodes.size() - 1; i >= 0; --i)
				nodes[i]->backward();

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->momentumSGD(0.1, 0.01);
			}

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->clear();
			}

			cout << n << "\t\t" << end->averageerror << "\t" << a->value << "\t" << endl;
		}


	}

	void test2() {
		inputnode * x1 = new inputnode(start);
		inputnode * x2 = new inputnode(start);
		x1->name = "x1";
		x2->name = "x2";
		doublevector x1data, x2data, y1data, y2data;
		x1data.push_back(1);
		x1data.push_back(2);
		x1data.push_back(3);
		//x1data.push_back(4);

		x2data.push_back(2);
		x2data.push_back(3);
		x2data.push_back(5);
		//x2data.push_back(7);

		y1data.push_back(5);
		y1data.push_back(7);
		y1data.push_back(9);
		//y1data.push_back(9);

		y2data.push_back(3);
		y2data.push_back(2);
		y2data.push_back(6);
		//y2data.push_back(6);

		inputnode * y1 = new inputnode(start);
		inputnode * y2 = new inputnode(start);
		start->data[x1->id] = x1data;
		start->data[x2->id] = x2data;
		start->data[y1->id] = y1data;
		start->data[y2->id] = y2data;
		LinearTransform * LT = new LinearTransform(2, 2);
		packnode * pn = new packnode;
		link(x1, pn);
		link(x2, pn);
		LT->test();
		link(pn, LT);
		unpacknode * up = new unpacknode;
		link(LT, up);

		addnode * t1 = new addnode(up);
		addnode * t2 = new addnode(up);
		t1->name = "t1";
		t2->name = "t2";

		subtractnode * e1 = sub(t1, y1);
		subtractnode * e2 = sub(t2, y2);

		endnode * end = new endnode(add(multi(e1, e1), multi(e2, e2)));

		topologysort();

		int n = 10000;
		while (n--) {
			for (int i = 0; i < nodes.size(); ++i)
				nodes[i]->forward();

			for (int i = nodes.size() - 1; i >= 0; --i)
				nodes[i]->backward();

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->momentumSGD(0.9, 0.01);
			}

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->clear();
			}

			cout << n << "\t" << end->averageerror << endl;
		}

		LT->explify();

		return;
	}
	void test() {
		//		inputnode * x = new inputnode(start);
		inputnode * y = new inputnode(start);
		inputnode * x = new inputnode();
		x->name = "x";
		y->name = "y";
		variablenode * a = new variablenode();
		a->init();
		a->name = "a";
		//multiplynode * multi = new multiplynode(x, a);
		variablenode * b = new variablenode();
		b->init();
		b->name = "b";
		//addnode * add = new addnode(b, multi);
		//subtractnode * sub = new subtractnode(add, y);
		add(b, multi(x, a));
		complexnode * ca = new complexnode(x);
		link(start, ca);
		abstractnode * error = abs(sub(y, ca));
		endnode * end = new endnode(error);
		doublevector xdata, ydata;
		xdata.push_back(1);
		xdata.push_back(2);
		xdata.push_back(3);
		xdata.push_back(4);
		ydata.push_back(2);
		ydata.push_back(3);
		ydata.push_back(5);
		ydata.push_back(7);
		start->data[ca->id] = xdata;
		start->data[y->id] = ydata;
		topologysort();
		int n = 3000;
		while (n--) {
			for (int i = 0; i < nodes.size(); ++i)
				nodes[i]->forward();

			for (int i = nodes.size() - 1; i >= 0; --i)
				nodes[i]->backward();

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->momentumSGD(0.1, 0.01);
			}

			for (int i = 0; i < nodes.size(); ++i) {
				nodes[i]->clear();
			}

			cout << n << "\t\t" << end->averageerror << "\t" << a->value << "\t" << b->value << endl;
		}

	}
};

int main()
{

	NN n;
	n.test2();
	return 0;
}

node::node(int _type, node * p) :type(_type) {
	switch (type) {
	case START:
		break;
	case ADD:
		middle = "+";
		break;
	case SUB:
		middle = "-";
		break;
	case MULTI:
		middle = "*";
		break;
	case ABS:
		left = "|";
		right = "|";
		break;
	}
	name = left + p->name + right;
	id = currentid++;
	NN::link(p, this);
}

node::node(int _type, node * p, node * q) :type(_type) {
	switch (type) {
	case START:
		break;
	case ADD:
		middle = "+";
		break;
	case SUB:
		middle = "-";
		break;
	case MULTI:
		middle = "*";
		break;
	case ABS:
		left = "|";
		right = "|";
		break;
	}
	name = "(" + p->name + ")" + middle + "(" + q->name + ")";
	id = currentid++;
	NN::link(p, this);
	NN::link(q, this);
}

node::node(int _type, vector<node *> _v) :type(_type) {
	switch (type) {
	case START:
		break;
	case ADD:
		middle = "+";
		break;
	case SUB:
		middle = "-";
		break;
	case MULTI:
		middle = "*";
		break;
	case ABS:
		left = "|";
		right = "|";
		break;
	}
	for (node * p : _v) {
		name = name + "(" + p->name + ")" + middle;
		NN::link(p, this);
	}
	name = name.substr(0, name.size() - 1);
}

inline addnode * NN::add(node * p, node * q) {
	return new addnode(p, q);
}

inline multiplynode * NN::multi(node * p, node * q) {
	return new multiplynode(p, q);
}

inline abstractnode * NN::abs(node * p) {
	return new abstractnode(p);
}

inline subtractnode * NN::sub(node * p) {
	return new subtractnode(p);
}

inline subtractnode * NN::sub(node * p, node * q) {
	return new subtractnode(p, q);
}

void NN::topologysort(node * p, vector<node *> & v) {
	v.push_back(p);

	for (int i = 0; i < p->succ_nodes.size(); ++i) {
		if (find(v.begin(), v.end(), p->succ_nodes[i]->succ) != v.end()) {
			continue;
		}
		node * q = p->succ_nodes[i]->succ;
		bool flag = true;
		for (int j = 0; j < q->pre_nodes.size(); ++j) {
			if (find(v.begin(), v.end(), q->pre_nodes[j]->pre) == v.end())
				flag = false;
		}
		if (flag) {
			topologysort(q, v);
		}
	}
}

void NN::findallnodes(node * p, vector<node*>& v)
{
	v.push_back(p);
	for (int i = 0; i < p->succ_nodes.size(); ++i) {
		if (find(v.begin(), v.end(), p->succ_nodes[i]->succ) != v.end()) {
			continue;
		}
		node * q = p->succ_nodes[i]->succ;
		findallnodes(q, v);
	}
	for (int i = 0; i < p->pre_nodes.size(); ++i) {
		if (find(v.begin(), v.end(), p->pre_nodes[i]->pre) != v.end()) {
			continue;
		}
		node * q = p->pre_nodes[i]->pre;
		findallnodes(q, v);
	}
}

vector<node*> NN::topologysort(node * p)
{
	vector<node *> allnodes, ans;
	findallnodes(p, allnodes);
	startnode temp;
	for (node * p : allnodes) {
		connector * c = new connector();
		c->succ = p;
		temp.succ_nodes.push_back(c);
	}
	topologysort(&temp, ans);
	ans.erase(ans.begin());
	return ans;
}

inline NN::NN() {
	start = new startnode();
}

inline void NN::link(node * p, node * q) {
	connector * c = new connector();
	c->pre = p;
	c->succ = q;
	p->succ_nodes.push_back(c);
	q->pre_nodes.push_back(c);
}

complexnode::complexnode(node * s) : node(COMPLEX)
{
	v = NN::topologysort(s);
	if (v.front() != s)
		throw(string("topologysort error"));
}

inline void complexnode::build(node * s) {
	v = NN::topologysort(s);
	if (v.front() != s)
		throw(string("topologysort error"));
}

inline complexnode * complexnode::copy() {
	vector<node *> temp;
	for (int i = 0; i < v.size(); ++i) {
		temp.push_back(v[i]->copy());
	}
	for (int i = 0; i < v.size(); ++i) {
		for (int j = 0; j < v[i]->succ_nodes.size(); ++j) {
			for (int k = 0; k < v.size(); ++k) {
				if (v[i]->succ_nodes[j]->succ == v[k]) {
					NN::link(temp[i], temp[k]);
					break;
				}
			}
		}
	}
	complexnode * ans = new complexnode(temp[0]);
	ans->name = this->name;
	return ans;
}

inline LinearTransform::LinearTransform(int _n, int _m) :n(_n), m(_m) {
	stringstream ss;
	ss << "Linear Transform " << n << string("*") << m;
	name = ss.str();
	ss.str("");
	packnode * out = new packnode;
	unpacknode * in = new unpacknode;
	x = new addnode[n];
	for (int i = 0; i < n; ++i) {
		ss << "x" << i + 1;
		x[i].name = ss.str();
		ss.str("");
		NN::link(in, x + i);
	}
	y = new addnode[m];
	for (int j = 0; j < m; ++j) {
		ss << "y" << j + 1;
		y[j].name = ss.str();
		ss.str("");
		NN::link(y + j, out);
	}
	a = new variablenode[n*m];
	b = new variablenode[n*m];
	for (int j = 0; j < m; ++j) {
		for (int i = 0; i < n; ++i) {
			ss << "a[" << i + 1 << "][" << j + 1 << "]";
			a[j*n + i].name = ss.str();
			ss.str("");
			NN::link(NN::multi(a + j*n + i, x + i), y + j);
		}
		ss << "b[" << j + 1 << "]";
		b[j].name = ss.str();
		ss.str("");
		NN::link(b + j, y + j);
	}
	build(in);
}

inline parallelfunction::parallelfunction(int _n, node * _t_n) : n(_n), template_node(_t_n) {
	template_node->clear();
	if (template_node->pre_nodes.size() != 1 || template_node->succ_nodes.size() != 1) {
		throw(string("not one to one"));
	}
	for (int i = 0; i < n; ++i)
		v.push_back(template_node->copy());
	unpacknode * up = new unpacknode;
	packnode * pn = new packnode;
	for (int i = 0; i < n; ++i) {
		NN::link(up, v[i]);
		NN::link(v[i], pn);
	}
	build(up);
	name = template_node->name;
}
