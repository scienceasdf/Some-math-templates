#ifndef SWH_SQueue_H_INCLUDED
#define SWH_SQueue_H_INCLUDED

template <class Type> class SQueue{
public:
    //empty SQueue
    SQueue():head(0),tail(0) {}
    //copy control to manage pointers to SQueueItems in the SQueue
    SQueue(const SQueue &Q):head(0),tail(0) {copy_elems(Q); }
    SQueue& operator=(const SQueue&);
    ~SQueue() {destroy(); }
    //return element from head of SQueue
    //unchecked operation:front on an empty SQueue is undefined
    Type& front()
    const Type &fornt() const {return head->item; }
    void push(const Type &);        //add element from head of SQueue
    void pop();                     //remove element from head of SQueue
    bool empty () const {           //true if no elements in the SQueue
        return head==0;
    }
private:
    SQueueItem<Type> *head;          //pointer to first element in SQueue
    SQueueItem<Type> *tail;          //pointer to last element in SQueue
    //utility functions used by copy constructor, assignment, and destructor
    void destroy();                 //delete all the elements
    void copy_elems(const SQueue&);  //copy elements from parameter
}
template <class Type> void SQueue<Type>::pop()
{
    //pop is unchecked: Popping off an empty SQueue id undefined
    SQueueItem<Type>* p=head;        //keep pointer to head so we can delete it
    head=head->next;                //head now points to next element
    delete p;
}

template <class Type> void SQueue<Type>::destroy()
{
    while(!empty())
        pop();
}

template <class Type> void SQueue<Type>::push(const type &val)
{
    //allocate a new SQueueItem object
    SQueueItem<Type> *pt=new SQueueItem<Type>(val);
    //put item onto existing SQueue
    if(empty())
        head=tail=pt;           //the SQueue now has only one element
    else {
        tail->next=pt;          //add new element to end of the SQueue
        tail=pt;
    }
}

template <class,Type> void SQueue<Type>::copy_elems(const SQueue &orig)
{
    //copy elements from orig into this SQueue
    //loop stops when pt==0, which happens when we reach orig.tail
    for(SQueueItem<type> *pt=orig.head;pt;pt=pt->next)
        push(pt->item);         //copy the element
}

SQueue& operator=(const SQueue& orig)
{
    for(SQueueItem<type> *pt=orig.head;pt;pt=pt->next)
        push(pt->item);
    return *this;
}

#endif // SWH_SQueue_H_INCLUDED
