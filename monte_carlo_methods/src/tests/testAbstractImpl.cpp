//created by John Stanco on 5/23/18

//testing implementation of methods for abstract class

#include <iostream>

class Base{
public:
	Base(){}
	virtual void A() = 0;
	void B(){printf("B\n");}

	virtual Base* C() = 0;
	void D(Base& b){b.A();}
	virtual Base& E() = 0;

};


class Derived : public Base{
protected:
	int n;
public:
	Derived() : n(0){}
	Derived(int n) : n(n){}
	void A(){printf("%d\n", this->n);}
	Derived* C(){return new Derived(this->n);}
	void D(Derived &d){d.A();}
	Derived& E(){return *this->C();}
};


int main(){
	Derived d(4);
	d.A();
	d.B();
	d.D(d);
	return 1;
}