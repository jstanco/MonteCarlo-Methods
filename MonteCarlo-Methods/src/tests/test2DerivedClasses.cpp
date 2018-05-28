//created by John Stanco on 5/17/18

#include <iostream>

class Base{
public:
	virtual void A() = 0;
};


class Derived : public Base{
	int n;
public:
	Derived() : n(0){};
	Derived(int n) : n(n){}
	void A(){printf("Derived\n");}
};


class usesBase{
	Base *b;
public:
	usesBase(Base *b) : b(b){}
	void B(){b->A();}
};


class usesDerived : public Derived{
	Derived *d;
public:
	usesDerived(Derived *d) : d(d){}
	void B(){d->A();}
};



int main(){
	Derived *d = new Derived(10);
	usesDerived uD(d);
	usesBase uB(d);
	uB.B();
	uD.B();

	return 1;
}