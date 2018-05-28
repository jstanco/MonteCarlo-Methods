//testing implementatino of derived class
#include <iostream>


class base{
public:
	virtual void A(){std::cout << "base" << std::endl;};
};


class derived : public base{
public:
	//derived(){}
	void A(){std::cout << "derived" << std::endl;}
};


class Foo{
public:
	base b;
	void B(){b.A();}
};


class Bar{
public:
	derived der;
	void B(){der.A();}
};


int main(){
	Foo foo;
	foo.B();
	Bar bar;
	bar.B();
}