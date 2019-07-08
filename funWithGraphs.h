#pragma once

#include <math.h>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class Node
{
protected:

    vector<shared_ptr<Node>> myArguments;

    bool        myProcessed = false;
    unsigned    myOrder = 0;
    double      myResult;

public:

    virtual ~Node() {}

    //  visitFunc: 
    //  a function of Node& that conducts a particular form of visit
    //  templated so we can pass a lambda
    template <class V>
    void postorder(V& visitFunc)
    {
        //  Already processed -> do nothing
        if (myProcessed == false)
        {
            //  Process ancestors first
            for (auto argument : myArguments)
                argument->postorder(visitFunc);

            //  Visit the node
            visitFunc(*this);
            
            //  Mark as processed
            myProcessed = true;
        }
    }
    
    //  visits
    
    virtual void evaluate() = 0;

    virtual void logInstruction() = 0;

    void setOrder(unsigned order)
    {
        myOrder = order;
    }

    //  Access result

    unsigned order()
    {
        return myOrder;
    }

    double result()
    {
        return myResult;
    }

    //  Reset processed flags
    void resetProcessed()
    {
        for (auto argument : myArguments) argument->resetProcessed();
        myProcessed = false;
    }
};

class PlusNode : public Node
{
public:

    PlusNode(shared_ptr<Node> lhs, shared_ptr<Node> rhs)
    {
        myArguments.resize(2);
        myArguments[0] = lhs;
        myArguments[1] = rhs;
    }

    void evaluate() override
    {
        myResult = myArguments[0]->result() + myArguments[1]->result();
    }

    void logInstruction() override
    {
        cout << "y" << order() 
            << " = y" << myArguments[0]->order() 
            << " + y" << myArguments[1]->order()
            << endl;
    }
};

class TimesNode : public Node
{
public:

    TimesNode(shared_ptr<Node> lhs, shared_ptr<Node> rhs)
    {
        myArguments.resize(2);
        myArguments[0] = lhs;
        myArguments[1] = rhs;
    }

    void evaluate() override
    {
        myResult = myArguments[0]->result() * myArguments[1]->result();
    }

    void logInstruction() override
    {
        cout << "y" << order()
            << " = y" << myArguments[0]->order()
            << " * y" << myArguments[1]->order()
            << endl;
    }
};

class LogNode : public Node
{
public:

    LogNode(shared_ptr<Node> arg)
    {
        myArguments.resize(1);
        myArguments[0] = arg;
    }

    void evaluate() override
    {
        myResult = log(myArguments[0]->result());
    }

    void logInstruction() override
    {
        cout << "y" << order() << " = log(" 
            << "y" << myArguments[0]->order() << ")" << endl;
    }
};

class Leaf: public Node
{

    double myValue;

public:

    Leaf(double val)
        : myValue(val) {}

    double getVal()
    {
        return myValue;
    }

    void setVal(double val)
    {
        myValue = val;
    }

    void evaluate() override
    {
        myResult = myValue;
    }

    void logInstruction() override
    {
        cout << "y" << order() << " = " << myValue << endl;
    }
};

class Number
{
    shared_ptr<Node> myNode;

public:

    Number(double val)
        : myNode(new Leaf(val)) {}

    Number(shared_ptr<Node> node)
        : myNode(node) {}

    shared_ptr<Node> node()
    {
        return myNode;
    }

    void setVal(double val)
    {
        //  Cast to leaf, only leaves can be changed
        dynamic_pointer_cast<Leaf>(myNode)->setVal(val);
    }

    double getVal()
    {
        //  Same comment here, only leaves can be read
        return dynamic_pointer_cast<Leaf>(myNode)->getVal();
    }

    double evaluate()
    {
        myNode->resetProcessed();
        auto visit = [](Node& n) {n.evaluate(); };
        myNode->postorder(visit);
        return myNode->result();
    }

    void setOrder()
    {
        myNode->resetProcessed();
        auto visit = [order = 0](Node& n) mutable {n.setOrder(++order); };
        myNode->postorder(visit);
    }

    void logResults()
    {
        myNode->resetProcessed();
        auto visit = [] (Node& n) 
        {
            cout << "Processed node " 
                << n.order() << " result = " 
                << n.result() << endl;
        };
        myNode->postorder(visit);
    }

    void logProgram()
    {
        myNode->resetProcessed();
        auto visit = [](Node& n) {n.logInstruction(); };
        myNode->postorder(visit);
    }
};

shared_ptr<Node> operator+(Number lhs, Number rhs)
{
    return shared_ptr<Node>(new PlusNode(lhs.node(), rhs.node()));
}

shared_ptr<Node> operator*(Number lhs, Number rhs)
{
    return shared_ptr<Node>(new TimesNode(lhs.node(), rhs.node()));
}

shared_ptr<Node> log(Number arg)
{
    return shared_ptr<Node>(new LogNode(arg.node()));
}

/*

inline double f(double x[5])
{
    double y1 = x[2] * (5.0 * x[0] + x[1]);
    double y2 = log(y1);
    double y = (y1 + x[3] * y2) * (y1 + y2);
    return y;
}

*/

template <class T>
inline T f(T x[5])
{
    auto y1 = x[2] * (5.0 * x[0] + x[1]);
    auto y2 = log(y1);
    auto y = (y1 + x[3] * y2) * (y1 + y2);
    return y;
}

void compute()
{
    //  Set inputs
    Number x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
       
    //  Build the dag
    Number y = f(x);

    //  Set order on the dag
    y.setOrder();

    //  Evaluate on the dag
    cout << y.evaluate() << endl;   // 797.751

    //  Log all results
    y.logResults();

    //  Log program 
    y.logProgram();

    //  Change x0 on the dag
    x[0].setVal(2.5);

    //  Evaluate the dag again
    cout << y.evaluate() << endl;   // 2769.76

    //  Log results again
    y.logResults();

    //  Log program 
    y.logProgram();  
}