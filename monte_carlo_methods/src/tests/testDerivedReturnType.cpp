//testing covariant return types

class Base {
public:
    virtual ~Base() {}
    virtual Base* clone() const;
};

class Derived : public Base {
public:
    virtual Derived* clone() const {
        return new Derived(*this);
    }
};

int main(){
	return 1;
}