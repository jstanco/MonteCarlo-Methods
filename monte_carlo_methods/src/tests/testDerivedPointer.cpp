//testing use of pointers to derived class
#include <iostream>


class base{
public:
	virtual void A(base*) = 0;
};




template<class T>
class derivedTemplate : public base{
protected:
	T data;
public:
	derivedTemplate() : data(T()){}
	derivedTemplate(const T &data) : data(data){}
	void A(base*) {std::cout << "derivedTemplate" << std::endl;}
};


template<class T>
class derived : public derivedTemplate<T>{
protected:
	int n;
public:
	derived() : n(0) {}
	derived(int n) : n(n) {}
	virtual void A(base*){std::cout << "derived" << std::endl;}
};


class Foo{
public:
	base* b;
	void B(){}
};


class Bar{
public:
	base* b;
	Bar(base* b) : b(b) {}
	void B(){}
};


int main(){
	derivedTemplate<int> d1(4);
	derivedTemplate<int> d2(5);
	derived<int> d(4);
	d2.A(&d1);
	return 1;
};