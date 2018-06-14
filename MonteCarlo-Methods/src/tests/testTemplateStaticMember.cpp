//created by John Stanco 6.1.18

#include <iostream>

template<char... Char>
class A{
	static int n;
};



int main(){
	A<'test'> a;
	
	return 1;
}