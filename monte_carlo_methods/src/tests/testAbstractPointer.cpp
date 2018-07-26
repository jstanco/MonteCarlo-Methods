//testing use of pointers to abstract class
#include <iostream>


class base{
public:
	virtual void A() = 0;
};


class Foo{
public:
	base* b;
	void B(){b->A();}
};


class Bar{
public:
	base* b;
	Bar(base* b) : b(b) {}
	void B(){b->A();}
};


int main(){
	Foo foo;
	foo.B();
	base* b;
	Bar bar(b);
	bar.B();


	return 1;
};