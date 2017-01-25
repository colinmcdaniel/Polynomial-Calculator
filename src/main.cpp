/*
 Colin McDaniel
 
 January 16, 2013
 
 This program prompts the user to input two polynomials, displays the polynomials, and then proceeds to perform various manipulations on these functions.
 */


#include<iostream>
#include<vector>
#include<cmath>

using namespace std;


class Polynomial {
    
public:
    
    // Constructors
    Polynomial();
    Polynomial(vector<int> coeffs);
    
    // Accessors
    int Degree() const;
    int Coefficient(int k) const;
    double evaluateAt(double x) const;
    void print() const;
    
    // Mutators
    void constantMultiply(int x);
    void Transform();
    
    // Member operators
    Polynomial operator++(int unused);
    Polynomial& operator++();
    Polynomial operator--(int unused);
    Polynomial& operator--();
    Polynomial& operator+=(Polynomial poly2);
    Polynomial& operator-=(Polynomial poly2);
    Polynomial& operator*=(Polynomial poly2);
    
private:
    vector<int> coefficient;
    
};


class Rational {
    
public:
    
    // Constructors
    Rational();
    Rational(Polynomial p);
    Rational(Polynomial pN, Polynomial pD);
    
    // Accessors
    double evaluateAt(double x);
    Polynomial getTop() const;
    Polynomial getBottom() const;
    void printRational();
    
    // Member operators
    Rational operator++(int unused);
    Rational& operator++();
    Rational operator--(int unused);
    Rational& operator--();
    Rational& operator+=(Rational rat2);
    Rational& operator-=(Rational rat2);
    Rational& operator*=(Rational rat2);
    Rational& operator/=(Rational rat2);
    
private:
    Polynomial poly_top;
    Polynomial poly_bottom;
    
};


// Nonmember functions that perform operations on polynomials
Polynomial Add(Polynomial poly1, Polynomial poly2);
Polynomial Subtract(Polynomial poly1, Polynomial poly2);
Polynomial Multiply(Polynomial poly1, Polynomial poly2);

// Nonmember operators for Polynomials
Polynomial operator+ (const Polynomial &poly1, const Polynomial &poly2);
Polynomial operator- (const Polynomial &poly1, const Polynomial &poly2);
Polynomial operator* (const Polynomial &poly1, const Polynomial &poly2);

// Nonmember operators for Rationals
Rational operator+ (const Rational &rat1, const Rational &rat2);
Rational operator- (const Rational &rat1, const Rational &rat2);
Rational operator* (const Rational &rat1, const Rational &rat2);
Rational operator/ (const Rational &rat1, const Rational &rat2);

// Boolean operators
bool operator==(Polynomial poly1, Polynomial poly2);
bool operator!=(Polynomial poly1, Polynomial poly2);
bool operator<=(Polynomial poly1, Polynomial poly2);
bool operator>=(Polynomial poly1, Polynomial poly2);
bool operator<(Polynomial poly1, Polynomial poly2);
bool operator>(Polynomial poly1, Polynomial poly2);


int main() {
    
    // Prints introduction to screen
    cout << "Welcome! Please input the coefficients of the first polynomial, p.\nWhen you are finished, enter -123456789.\n";
    
    // Declares vector v1 for storage
    int coeffs_input1 = 0;
    vector<int> v1;
    
    // Gathers first polynomial data from user and stores in vector v1
    while(coeffs_input1 != -123456789){
        cin >> coeffs_input1;
        
        if(coeffs_input1 != -123456789)
            v1.push_back(coeffs_input1);
    }
    
    // Create Polynomial p from v1
    Polynomial p(v1);
    
    // Prints Polynomial p to screen
    cout << endl << "Your first polynomial is p(x) = ";
    p.print();
    
    // Prompts user to input second polynomial
    cout << endl << endl << "Please input the coefficients of the second polynomial, q.\n";
    
    // Declares vector v2 for storage
    int coeffs_input2 = 0;
    vector<int> v2;
    
    // Gathers second polynomial data from user and stores in vector v2
    while(coeffs_input2 != -123456789){
        cin >> coeffs_input2;
        
        if(coeffs_input2 != -123456789)
            v2.push_back(coeffs_input2);
    }
    
    // Create Polynomial q from v2
    Polynomial q(v2);
    
    // Prints Polynomial q to screen
    cout << endl << "Your second polynomial is q(x) = ";
    q.print();
    
    // Outputs p + q
    cout << endl << endl << "p(x)+q(x) = ";
    (p + q).print();
    
    // Outputs p - q
    cout << endl << endl << "p(x)-q(x) = ";
    (p - q).print();
    
    // Outputs p * q
    cout << endl << endl << "p(x)*q(x) = ";
    (p * q).print();
    
    // Outputs the Rational (p / q)
    cout << endl << endl << "p(x)/q(x) = ";
    Rational r_numPdenQ(p, q);
    r_numPdenQ.printRational();
    
    // Outputs the sum of the Rationals (p / q) + (p * q)
    cout << endl << endl << "p(x)/q(x)+p(x)*q(x) = ";
    Rational r_numPQ(p * q);
    (r_numPdenQ + r_numPQ).printRational();
    
    // Outputs the polynomial p incremented by 1
    cout << endl << endl << "p(x)+1 = ";
    (++p).print();
    
    // Outputs the polynomial p incremented by 1 again
    cout << endl << endl << "p(x)+2 = ";
    (++p).print();
    
    // Outputs the product of the Rationals (p / q) * ((1+x^2-3x^4) / 1)
    cout << endl << endl << "(p(x)/q(x)) * (1+x^2-3x^4) = ";
    
    int degree0 = 1;
    int degree1 = 0;
    int degree2 = 1;
    int degree3 = 0;
    int degree4 = -3;
    
    vector<int> v_given (5);
    v_given[0] = degree0;
    v_given[1] = degree1;
    v_given[2] = degree2;
    v_given[3] = degree3;
    v_given[4] = degree4;
    
    Polynomial p_given(v_given);
    Rational r_given(p_given);
    (r_numPdenQ * r_given).printRational();
    
    // Decrements Polynomial p to its original state for testing equality to q
    p--;
    p--;
    
    // Outputs whether p equals q
    cout << endl << endl << "Does p(x) equal q(x)? ";
    if(p == q)
        cout << "Yes.";
    else
        cout << "No.";
    
    // Outputs whether p is less than or equal to q
    cout << endl << endl << "Is p(x) <= q(x)? ";
    if(p <= q)
        cout << "Yes.";
    else
        cout << "No.";
    
    cout << endl << endl;
    
    return 0;
}


/** Polynomial default constructor, initializes polynomial to 0
 */
Polynomial::Polynomial() {
    coefficient.push_back(0);
}


/** Polynomial constructor
 @param coeffs Value gets set as private variable coefficient
 */
Polynomial::Polynomial(vector<int> coeffs) {
    coefficient = coeffs;
}


/** Polynomial member function definition
 @return Degree of polynomial
 */
int Polynomial::Degree() const {
    int degree = static_cast<int>(coefficient.size() - 1);
    
    return degree;
}


/** Polynomial member function definition, returns the coefficient of x^k of a given k
 @param k Degree of polynomial
 @return The coefficient of the given degree
 */
int Polynomial::Coefficient(int k) const {
    return coefficient[k];
}


/** Polynomial member function definition, returns the evaluation of the polynomial at a given x
 @param x Where polynomial is evaluated
 @return f(x)
 */
double Polynomial::evaluateAt(double x) const {
    double counter = 0;
    
    for(int i = 0; i < coefficient.size(); i++)
        counter = coefficient[i] * pow(x, i) + counter;
    
    return counter;
}


/** Polynomial member function definition, prints the polynomial to screen
 */
void Polynomial::print() const {
    
    // Marks running total of coefficients in the for loop
    int sum_tally = 0;
    
    // Represents the sum of all the coefficients
    int sum_elements = 0;
    
    // Checks if given coefficient is the first nonzero coefficient in the polynomial
    int zero_marker = 0;
    
    // Tallies when running total is not zero in the for loop
    int tally = 0;
    
    // Sums all coefficients
    for(int i = 0; i < coefficient.size(); i++)
        sum_elements += coefficient[i];
    
    // Prints polynomial to screen
    for(int i = 0; i < coefficient.size(); i++) {
        if(coefficient[i] == 0){
            if(sum_elements == 0 && zero_marker == 0) {
                cout << 0;
                
                // Marks that 0 has been printed to prevent this from happening more than once during later iterations
                zero_marker++;
            }
            else
                cout << "";
        }
        else{
            if(i == 0)
                cout << coefficient[i];
            else{
                if(tally != 0 && coefficient[i] > 0)
                    cout << "+";
                
                if(i == 1) {
                    if(coefficient[i] < -1)
                        cout << coefficient[i] << "x";
                    else if(coefficient[i] == -1)
                        cout << "-x";
                    else if(coefficient[i] == 1)
                        cout << "x";
                    else
                        cout << coefficient[i] << "x";
                }
                else{
                    if(coefficient[i] < -1)
                        cout << coefficient[i] << "x^" << i;
                    else if(coefficient[i] == -1)
                        cout << "-x^" << i;
                    else if(coefficient[i] == 1)
                        cout << "x^" << i;
                    else
                        cout << coefficient[i] << "x^" << i;
                }
            }
        }
        
        // Mark running total of coefficients
        sum_tally += coefficient[i];
        
        // Marks that the running total is not zero anymore (for plus sign)
        if(sum_tally != 0)
            tally++;
    }
    
    return;
}


/** Polynomial member function definition, multiplies polynomial by a given integer
 @param x Constant multiplied to polynomial
 */
void Polynomial::constantMultiply(int x) {
    for(int i = 0; i < coefficient.size(); i++)
        coefficient[i] *= x;
    
    return;
}


/** Polynomial member function definition, modifies polynomial as described (takes derivative)
 */
void Polynomial::Transform() {
    for(int i = 0; i < coefficient.size(); i++) {
        if(coefficient.size() == 1)
            coefficient[i] = 0;
        else if(coefficient.size() > 1 && i == coefficient.size() - 1)
            coefficient.pop_back();
        else
            coefficient[i] = (i + 1) * coefficient[i + 1];
    }
    
    return;
}


/** Polynomial postfix ++ member operator definition; copies Polynomial object, increments the original by one, and returns the copy
 @param unused Distinguishes operator as postfix
 @return Copy of Polynomial object
 */
Polynomial Polynomial::operator++(int unused) {
    Polynomial clone(coefficient);
    
    coefficient[0]++;
    
    return clone;
}


/** Polynomial prefix ++ member operator definition, increments Polynomial by 1
 @return Modified Polynomial object
 */
Polynomial& Polynomial::operator++() {
    coefficient[0]++;
    
    return *this;
}


/** Polynomial postfix -- member operator definition; creates copy of Polynomial object, decrements original by 1, and returns copy
 @param unused Distinguishes operator as postfix
 @return original Polynomial object
 */
Polynomial Polynomial::operator--(int unused) {
    Polynomial clone(coefficient);
    
    coefficient[0]--;
    
    return clone;
}


/** Polynomial prefix -- member operator definition; decrements Polynomial by 1
 @return Modified Polynomial object
 */
Polynomial& Polynomial::operator--() {
    coefficient[0]--;
    
    return *this;
}


/** Polynomial += member operator definition, adds the new polynomial to the current Polynomial and assigns the sum to that object
 @param poly2 Added to the current Polynomial object
 @return Modified Polynomial object
 */
Polynomial& Polynomial::operator+=(Polynomial poly2) {
    if((coefficient.size() - 1) < poly2.Degree()) {
        for(int i = 0; i < coefficient.size(); i++)
            coefficient [i] += poly2.Coefficient(i);
        
        for(int i = 0; i < ((poly2.Degree() + 1) - coefficient.size()); i++)
            coefficient.push_back(poly2.Coefficient(static_cast<int>(i + coefficient.size())));
    }
    else {
        for(int i = 0; i <= poly2.Degree(); i++)
            coefficient [i] += poly2.Coefficient(i);
    }
    
    return *this;
}


/** Polynomial -= member operator definition, subtracts the new polynomial from the current Polynomial and assigns the difference to that object
 @param poly2 Subtracted from the current Polynomial object
 @return Modified Polynomial object
 */
Polynomial& Polynomial::operator-=(Polynomial poly2) {
    
    // poly2 becomes negative
    poly2.constantMultiply(-1);
    
    // Follow += operator code
    if((coefficient.size() - 1) < poly2.Degree()) {
        for(int i = 0; i < coefficient.size(); i++)
            coefficient [i] += poly2.Coefficient(i);
        
        for(int i = 0; i < ((poly2.Degree() + 1) - coefficient.size()); i++)
            coefficient.push_back(poly2.Coefficient(static_cast<int>(i + coefficient.size())));
    }
    else {
        for(int i = 0; i <= poly2.Degree(); i++)
            coefficient [i] += poly2.Coefficient(i);
    }
    
    return *this;
}


/** Polynomial *= member operator definition, multiplies the new polynomial from the current Polynomial and assigns the product to that object
 @param poly 2 Multiplied to the current Polynomial object
 @return Modified Polynomial object
 */
Polynomial& Polynomial::operator*=(Polynomial poly2) {
    
    // Declare clone vector with appropriate number of elements to hold product
    vector<int> clone(coefficient.size() + poly2.Degree(), 0);
    
    // Perform multiplication and store in clone
    for(int i = 0; i < coefficient.size(); i++) {
        for(int j = 0; j < (poly2.Degree() + 1); j++)
            clone[i+j] += coefficient[i] * poly2.Coefficient(j);
    }
    
    // Give private variable coefficient the values of clone
    coefficient = clone;
    
    return *this;
}


/** Rational default constructor, initializes rational function to 0/1 by default
 */
Rational::Rational() {
    vector<int> vec_top(1, 0);
    vector<int> vec_bottom(1, 1);
    
    Polynomial clone1(vec_top);
    Polynomial clone2(vec_bottom);
    
    poly_top = clone1;
    poly_bottom = clone2;
}


/**
 Rational constructor, sets numerator private variable equal to p and denominator private variable equal to 1
 @param p Numerator
 */
Rational::Rational(Polynomial p) {
    poly_top = p;
    
    vector<int> vec_bottom(1, 1);
    Polynomial clone(vec_bottom);
    
    poly_bottom = clone;
}


/** Rational constructor, sets private variable for numerator equal to pN and private variable for denominator equal to pD
 @param pN Numerator
 @param pD Denominator
 */
Rational::Rational(Polynomial pN, Polynomial pD) {
    poly_top = pN;
    poly_bottom = pD;
}


/** Rational member function definition, evaluates the Rational at a given x
 @param x Where rational should be evaluated
 @return f(x)
 */
double Rational::evaluateAt(double x) {
    return (poly_top.evaluateAt(x) / poly_bottom.evaluateAt(x));
}


/** Rational member function definition
 @return Numerator
 */
Polynomial Rational::getTop() const {
    return poly_top;
}


/** Rational member function definition
 @return Denominator
 */
Polynomial Rational::getBottom() const {
    return poly_bottom;
}


/** Rational member function definition, prints the Rational to screen
 */
void Rational::printRational() {
    poly_top.print();
    cout << " / ";
    poly_bottom.print();
    
    return;
}


/** Rational Postfix ++ member operator definition; copies Rational, increments original Rational by 1, and returns copy
 @param unused Distinguishes as postfix
 @return Copy of original Rational object
 */
Rational Rational::operator++(int unused) {
    Rational clone(poly_top, poly_bottom);
    
    poly_top += poly_bottom;
    
    return clone;
}


/** Rational prefix ++ member operator definition; increments Rational object by 1
 @return Modified Rational object
 */
Rational& Rational::operator++() {
    poly_top += poly_bottom;
    
    return *this;
}


/** Rational postfix -- member operator definition; copies Rational, decrements original by 1, and returns copy
 @param unused Distinguishes as postfix
 @return Copy of original Rational object
 */
Rational Rational::operator--(int unused) {
    Rational clone(poly_top, poly_bottom);
    
    poly_top -= poly_bottom;
    
    return clone;
}


/** Ratinoal prefix -- member operator definition; decrements Rational object by 1
 @return Modified Rational object
 */
Rational& Rational::operator--() {
    poly_top -= poly_bottom;
    
    return *this;
}


/** Rational += member operator definition, add rat2 to current Rational object
 @param rat2 Second rational being added
 @return Modified Rational object
 */
Rational& Rational::operator+=(Rational rat2) {
    Rational rat1(poly_top, poly_bottom);
    
    Rational s((rat1.getTop() * rat2.getBottom()) + (rat2.getTop() * rat1.getBottom()), rat1.getBottom() * rat2.getBottom());
    
    poly_top = s.getTop();
    poly_bottom = s.getBottom();
    
    return *this;
}


/** Rational -= member operator definition, subtracts rat2 from current Rational object
 @param rat2 Second rational being subtracted
 @return Modified Rational object
 */
Rational& Rational::operator-=(Rational rat2) {
    Rational rat1(poly_top, poly_bottom);
    
    Rational d((rat1.getTop() * rat2.getBottom()) - (rat2.getTop() * rat1.getBottom()), rat1.getBottom() * rat2.getBottom());
    
    poly_top = d.getTop();
    poly_bottom = d.getBottom();
    
    return *this;
}


/** Rational *= member operator definition, multiplies rat2 with current Rational object and stores product in that object
 @param rat2 Second rational being multiplied
 @return Modified Rational object
 */
Rational& Rational::operator*=(Rational rat2) {
    Rational rat1(poly_top, poly_bottom);
    
    Rational p(rat1.getTop() * rat2.getTop(), rat1.getBottom() * rat2.getBottom());
    
    poly_top = p.getTop();
    poly_bottom = p.getBottom();
    
    return *this;
}


/** Rational /= member operator definition, divides rat2 from current Rational object and stores quotient in that object
 @param rat2 Second rational being divided
 @return Modified Rational object
 */
Rational& Rational::operator/=(Rational rat2) {
    Rational rat1(poly_top, poly_bottom);
    
    Rational q(rat1.getTop() * rat2.getBottom(), rat2.getTop() * rat1.getBottom());
    
    poly_top = q.getTop();
    poly_bottom = q.getBottom();
    
    return *this;
}


/** Nonmember function definition, adds two polynomials together and returns a Polynomial object
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return poly1 + poly2
 */
Polynomial Add(Polynomial poly1, Polynomial poly2){
    vector<int> vector_sum;
    
    // Enter the sum up to the last element of the smaller vector to the new vector, then enter the remaining elements of the larger vector to the new vector
    if(poly1.Degree() < poly2.Degree()) {
        for(int i = 0; i <= poly1.Degree(); i++)
            vector_sum.push_back(poly1.Coefficient(i) + poly2.Coefficient(i));
        
        for(int i = 0; i < (poly2.Degree() - poly1.Degree()); i++)
            vector_sum.push_back(poly2.Coefficient(i + poly1.Degree() + 1));
    }
    else {
        for(int i = 0; i <= poly2.Degree(); i++)
            vector_sum.push_back(poly1.Coefficient(i) + poly2.Coefficient(i));
        
        for(int i = 0; i < (poly1.Degree() - poly2.Degree()); i++)
            vector_sum.push_back(poly1.Coefficient(i + poly2.Degree() + 1));
    }
    
    // Create Polynomial object using new vector
    Polynomial sum(vector_sum);
    
    return sum;
}


/** Nonmember function definition, subtracts two polynomials and returns a Polynomial object
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return poly1 - poly2
 */
Polynomial Subtract(Polynomial poly1, Polynomial poly2) {
    poly2.constantMultiply(-1);
    return Add(poly1, poly2);
}


/** Nonmember function definition, mulitplies two polynomials and returns a Polynomial object
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return poly1 * poly2
 */
Polynomial Multiply(Polynomial poly1, Polynomial poly2) {
    vector<int> vector_product(poly1.Degree() + poly2.Degree() + 1, 0);
    
    for(int i = 0; i < (poly1.Degree() + 1); i++) {
        for(int j = 0; j < (poly2.Degree() + 1); j++)
            vector_product[i+j] += poly1.Coefficient(i) * poly2.Coefficient(j);
    }
    
    Polynomial product(vector_product);
    
    return product;
}


/** Addition operator, adds two &Polynomials together and returns a polynomial
 @param &poly1 Left side polynomial
 @param &poly2 Right side polynomial
 @return &poly1 + &poly2 Sum Polynomial
 */
Polynomial operator+ (const Polynomial &poly1, const Polynomial &poly2) {
    return Add(poly1, poly2);
}


/** Subtraction operator, subtracts one &Polynomial from another and returns a polynomial
 @param &poly1 Left side polynommial
 @param &poly2 Right side polynomial
 @return &poly1 - &poly2 Difference Polynomial
 */
Polynomial operator- (const Polynomial &poly1, const Polynomial &poly2) {
    return Subtract(poly1, poly2);
}


/** Multiplication operator, multiplies two &Polynomials together and returns a polynomial
 @param &poly1 Left side polynomial
 @param &poly2 Right side polynomial
 @return &poly1 * &poly2 Product Polynomial
 */
Polynomial operator* (const Polynomial &poly1, const Polynomial &poly2) {
    return Multiply(poly1, poly2);
}


/** Addition operator, adds two &Rationals together and returns a Rational
 @param &rat1 Left side Rational
 @param &rat2 Right side Rational
 @return &rat1 + &rat2 Sum Rational
 */
Rational operator+ (const Rational &rat1, const Rational &rat2) {
    Rational s((rat1.getTop() * rat2.getBottom()) + (rat2.getTop() * rat1.getBottom()), rat1.getBottom() * rat2.getBottom());
    
    return s;
}


/** Subtraction operator, subtracts one &Rational from another and returns a Rational
 @param &rat1 Left side Rational
 @param &rat2 Right side Rational
 @return &rat1 - &rat2 Difference Rational
 */
Rational operator- (const Rational &rat1, const Rational &rat2) {
    Rational d((rat1.getTop() * rat2.getBottom()) - (rat2.getTop() * rat1.getBottom()), rat1.getBottom() * rat2.getBottom());
    
    return d;
}


/** Multiplication operator, multiplies two &Rationals together and returns a Rational
 @param &rat1 Left side Rational
 @param &rat2 Righ side Rational
 @return &rat1 * &rat2 Product Rational
 */
Rational operator* (const Rational &rat1, const Rational &rat2) {
    Rational p(rat1.getTop() * rat2.getTop(), rat1.getBottom() * rat2.getBottom());
    
    return p;
}


/** Division operator, divides one &Rational by another and returns a Rational
 @param &rat1 Left side Rational
 @param &rat2 Right side Rational
 @return &rat1 / &rat2 Quotient Rational
 */
Rational operator/ (const Rational &rat1, const Rational &rat2) {
    Rational q(rat1.getTop() * rat2.getBottom(), rat2.getTop() * rat1.getBottom());
    
    return q;
}


/** == boolean operator definition, determines if one polynomial is equal to another
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return True or False
 */
bool operator==(Polynomial poly1, Polynomial poly2) {
    bool val = 0;
    int counter = 0;
    
	// Checks if elements are equal up to the end of the smallest vector, then checks if all the remaining elements on the larger vector are zero
	if(poly1.Degree() > poly2.Degree()) {
		for(int i = 0; i <= poly2.Degree(); i++){
			if(poly1.Coefficient(i) == poly2.Coefficient(i))
				counter++;
		}
        
		if(counter == (poly2.Degree() + 1)){
			for(int j = 0; j < (poly1.Degree() - poly2.Degree()); j++){
				if(poly1.Coefficient(j + poly2.Degree() + 1) == 0)
					counter++;
			}
		}
        
		if(counter == (poly1.Degree() + 1))
            val = 1;
	}
    
	if(poly2.Degree() > poly1.Degree()) {
		for(int i = 0; i <= poly1.Degree(); i++){
			if(poly2.Coefficient(i) == poly1.Coefficient(i))
				counter++;
		}
        
		if(counter == (poly1.Degree() + 1)){
			for(int j = 0; j < (poly2.Degree() - poly1.Degree()); j++){
				if(poly2.Coefficient(j + poly1.Degree() + 1) == 0)
					counter++;
			}
		}
        
		if(counter == (poly2.Degree() + 1))
			val = 1;
	}
    
    if(poly1.Degree() == poly2.Degree()) {
        for(int i = 0; i <= poly1.Degree(); i++){
            if(poly1.Coefficient(i) == poly2.Coefficient(i))
                counter++;
        }
        
		if(counter == (poly1.Degree() + 1))
            val = 1;
    }
    
    return val;
}


/** != boolean operator definition, determines if one polynomial is not equal to another
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return True or False
 */
bool operator!=(Polynomial poly1, Polynomial poly2) {
    bool val = 1;
    
    if(poly1 == poly2)
        val = 0;
    
    return val;
}


/** <= boolean operator definiton, determines if one polynomial is less than or equal to another
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return True or False
 */
bool operator<=(Polynomial poly1, Polynomial poly2) {
    bool val;
    
    if(poly1.Degree() <= poly2.Degree())
        val = 1;
    else
        val = 0;
    
    return val;
}


/**
 >= boolean operator definition, determines if one polynomial is greater than or equal to another
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return True or False
 */
bool operator>=(Polynomial poly1, Polynomial poly2) {
    bool val;
    
    if(poly1.Degree() >= poly2.Degree())
        val = 1;
    else
        val = 0;
    
    return val;
}


/** < boolean operator definiton, determines if one polynomial is less than another
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return True or False
 */
bool operator<(Polynomial poly1, Polynomial poly2) {
    bool val;
    
    if(poly1.Degree() < poly2.Degree())
        val = 1;
    else
        val = 0;
    
    return val;
}


/** > boolean operator definition, determines if one polynomial is greater than another
 @param poly1 Left side polynomial
 @param poly2 Right side polynomial
 @return True or False
 */
bool operator>(Polynomial poly1, Polynomial poly2) {
    bool val;
    
    if(poly1.Degree() > poly2.Degree())
        val = 1;
    else
        val = 0;
    
    return val;
}