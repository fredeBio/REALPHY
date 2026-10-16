package util;

import java.util.ArrayList;
//IntervalTree class
//
// CONSTRUCTION: with no initializer
//
// ******************PUBLIC OPERATIONS*********************
// void insert( x )       --> Insert x
// void remove( x )       --> Remove x
// Comparable find( x )   --> Return item that matches x
// Comparable findMin( )  --> Return smallest item
// Comparable findMax( )  --> Return largest item
// boolean isEmpty( )     --> Return true if empty; else false
// void makeEmpty( )      --> Remove all items
// ******************ERRORS********************************
// Exceptions are thrown by insert and remove if warranted

/**
 * Implements an Interval-tree.
 * Note that all "matching" is based on the compareTo method.
 * @author Mark Allen Weiss
 */
public class IntervalTree <T extends Interval<T>> {
    /**
     * Construct the tree.
     */
    public IntervalTree( ) {
        root = nullNode;
    }
    
    public void search(T inter,ArrayList<T> al){
    	search(root,inter,al);
    }
    
    /**
     * Insert into the tree.
     * @param x the item to insert.
     * @throws DuplicateItemException if x is already present.
     */
    public void insert( T x ) {
    	
        root = insert( x, root,null );
//       	if(counter==15){
//       		parseTree(root);
//       		System.exit(-1);
//       	}
//			System.out.println("________________");

    }
    
    /**
     * Remove from the tree.
     * @param x the item to remove.
     * @throws ItemNotFoundException if x is not found.
     */
    public void remove( T x ) {
        deletedNode = nullNode;
        root = remove( x, root );
    }
    
    /**
     * Find the smallest item in the tree.
     * @return the smallest item or null if empty.
     */
    public T findMin( ) {
        if( isEmpty( ) )
            return null;
        
        IntervalNode<T> ptr = root;
        
        while( ptr.left != nullNode )
            ptr = ptr.left;
        
        return ptr.element;
    }
    
    /**
     * Find the largest item in the tree.
     * @return the largest item or null if empty.
     */
    public T findMax( ) {
        if( isEmpty( ) )
            return null;
        
        IntervalNode<T> ptr = root;
        
        while( ptr.right != nullNode )
            ptr = ptr.right;
        
        return ptr.element;
    }
	// Search for all intervals which contain "p", starting with the
	// node "n" and adding matching intervals to the list "result"
	public void search(IntervalNode<T> n, T p, ArrayList<T> result) {
	    //System.out.println(p.start +" "+n.element.start+" "+n.max+" "+n.element.overlapsWith(p));

		// Don't search nodes that don't exist
	    if (n == nullNode)
	        return; 

	    
	    // If p is to the right of the rightmost point of any interval
	    // in this node and all children, there won't be any matches.

	    if (p.start>n.max){
	    	return;
	    }
	    // Search left children
	    if (n.left != nullNode){

	    	search(n.left, p, result);
	    }
	    // Check this node
	    if (n.element.overlapsWith(p)){
	        result.add(n.element);
	    }

	    // If p is to the left of the start of this interval,
	    // then it can't be in any child to the right.
	    if (p.end<n.element.start)
	        return;

	    // Otherwise, search right children
	    if (n.right != nullNode)
	        search(n.right, p, result);

	}
	/**
     * Find the predecessor and successor of an element...only works if element does not exist in tree
     */

    public Flank<T> findFlank( T x ) {
        IntervalNode<T> current = root;
       
        //nullNode.element = x;
        Flank<T> f=new Flank<T>();
        for( ; ; ) {
            if(current.element!=null && x.compareTo( current.element ) < 0 ){
                f.next=current.element;
            	current = current.left;
            }
            else if(current.element!=null && x.compareTo( current.element ) > 0 ){
            	f.prev=current.element;
            	current = current.right;
            }
            else 
    
            	return f;
                
        }
    }
    
    public T maxTree(IntervalNode<T> n){
    	while(n.right!=nullNode)n=n.right;
    	return n.element;
    }
    public T minTree(IntervalNode<T> n){
    	while(n.left!=nullNode)n=n.left;
    	return n.element;
    }
    public T predecessor(T x){
    	IntervalNode<T> current = root;
    	for( ; ; ) {
    		if(current.element==null)return null;
    		if( x.compareTo( current.element ) < 0 )
    			current = current.left;
    		else if( x.compareTo( current.element ) > 0 )
    			current = current.right;
    		else if( current != nullNode ){

    			if(current.left!=nullNode){
    				return maxTree(current.left);
    			}else{
    				IntervalNode<T> p=current.dad;
    				while(p!=null&&current==p.left){
    					current=p;
    					p=p.dad;
    				}
    				if(p!=null)
    					return p.element;
    				else return null;
    			}
    		}else
    			return null;
    	}
    }
    public T successor(T x){
    	IntervalNode<T> current = root;
    	for( ; ; ) {
    		if(current.element==null)return null;
    		if( x.compareTo( current.element ) < 0 )
    			current = current.left;
    		else if( x.compareTo( current.element ) > 0 )
    			current = current.right;
    		else if( current != nullNode ){

    			if(current.right!=nullNode){
    				return minTree(current.right);
    			}else{
    				IntervalNode<T> p=current.dad;
    				while(p!=null&&current==p.right){
    					current=p;
    					p=p.dad;
    				}
    				if(p!=null)
    					return p.element;
    				else return null;
    			}
    		}else
    			return null;
    	}
    }
	/**
     * Find an item in the tree.
     * @param x the item to search for.
     * @return the matching item of null if not found.
     */

    public T find( T x ) {
        IntervalNode<T> current = root;
        //nullNode.element = x;
        
        for( ; ; ) {
        	if(current.element==null)return null;
            if( x.compareTo( current.element ) < 0 )
                current = current.left;
            else if( x.compareTo( current.element ) > 0 )
                current = current.right;
            else if( current != nullNode )
                return current.element;
            else
                return null;
        }
    }
    
    /**
     * Make the tree logically empty.
     */
    public void makeEmpty( ) {
        root = nullNode;
    }
    
    /**
     * Test if the tree is logically empty.
     * @return true if empty, false otherwise.
     */
    public boolean isEmpty( ) {
        return root == nullNode;
    }
    
    /*
     * Goes up to the top of the tree and adjust the right end of the interval for every node it passes
     */

    private void adjustMax(IntervalNode<T> t,int max){
    	if(t!=null){
    		if (t.max<max){
    			t.max=max;
    			adjustMax(t.dad,max);
    		}
    	}
    }
    /*
     * goes to the dad and checks if max == max if yes, max=currentMax or Max of left, then adjust max algorithm
     */

    private void adjustMaxDeletion(IntervalNode<T> t){
    	if(t!=null){
    			
    			t.max=max(t.element.getEnd(),t.left.max,t.right.max);
    			//System.out.println(t.element+" "+t.element.getEnd()+" "+t.left.max+" "+t.right.max+" "+t.right.element);
    			//adjustMaxDeletion(t.dad);

    	}
    }
    /**
     * Internal method to insert into a subtree.
     * @param x the item to insert.
     * @param t the node that roots the tree.
     * @return the new root.
     * @throws DuplicateItemException if x is already present.
     */
    private IntervalNode<T> insert( T x, IntervalNode<T> t,IntervalNode<T> father ) {
 
    	if( t == nullNode ){
            t = new IntervalNode<T>( x,father );
            adjustMax(t.dad,x.end);
        }
        else if( x.compareTo( t.element ) < 0 )
            t.left = insert( x, t.left,t );
        else if( x.compareTo( t.element ) > 0 )
            t.right = insert( x, t.right,t );
        else {t.element.append(x);
           // throw new DuplicateItemException( x.toString( ) );
        }
        t = skew( t );
        t = split( t );
        return t;
    }
    
    
    /**
     * Internal method to remove from a subtree.
     * @param x the item to remove.
     * @param t the node that roots the tree.
     * @return the new root.
     * @throws ItemNotFoundException if x is not found.
     * 
     * 
     */
    private IntervalNode<T> remove( Comparable<? super T> x, IntervalNode<T> t ) {
    	if( t != nullNode ) {
    		// Step 1: Search down the tree and set lastNode and deletedNode
    		lastNode = t;
    		if( x.compareTo( t.element ) < 0 ){
    			t.left = remove( x, t.left );
    			adjustMaxDeletion(t);
    		}
    		else {
    			deletedNode = t;
    			t.right = remove( x, t.right );
    			adjustMaxDeletion(t);
    		}

    		// Step 2: If at the bottom of the tree and
    		//         x is present, we remove it
    		if( t == lastNode ) {
    			if( deletedNode == nullNode || x.compareTo( deletedNode.element ) != 0 )
    				throw new ItemNotFoundException( x.toString( ) );
    			deletedNode.element = t.element;
    			
    			t = t.right;    			
    			
    		}

    		// Step 3: Otherwise, we are not at the bottom; rebalance
    		else
    			if( t.left.level < t.level - 1 || t.right.level < t.level - 1 ) {
    				if( t.right.level > --t.level )
    					t.right.level = t.level;

    				t = skew( t );
    				t.right = skew( t.right );
    				t.right.right = skew( t.right.right );
    				t = split( t );
    				t.right = split( t.right );
    			}
    		
    	}
    	
    	return t;
    }
    
    /**
     * Skew primitive for Interval-trees.
     * @param t the node that roots the tree.
     * @return the new root after the rotation.
     */
    private  IntervalNode<T> skew( IntervalNode<T> t ) {
        if( t.left.level == t.level )
            t = rotateWithLeftChild( t );
        return t;
    }
    
    /**
     * Split primitive for Interval-trees.
     * @param t the node that roots the tree.
     * @return the new root after the rotation.
     */
    private  IntervalNode<T> split( IntervalNode<T> t ) {
        if( t.right.right.level == t.level ) {
            t = rotateWithRightChild( t );
            t.level++;
        }
        return t;
    }
    
    /**
     * Rotate binary tree node with left child.
     */
    private  IntervalNode<T> rotateWithLeftChild( IntervalNode<T> k2 ) {
    	//System.out.println("rotateLeft");
    	//System.out.println("_____xxx");
    	//parseTree(root);
    	//System.out.println("_____xxx");
        IntervalNode<T> k1 = k2.left;
        k1.dad=k2.dad;
        k2.left = k1.right;
        k1.right = k2;
        k1.max=k2.max;
        //System.out.println(k1.left.left.max);
        //System.out.println(k1.left.right.max);
        //System.out.println(k1.left.element.end);

        int rightleft=k1.right.left==null?Integer.MIN_VALUE:k1.right.left.max;
        int rightright=k1.right.right==null?Integer.MIN_VALUE:k1.right.right.max;
        int rightelement=k1.right==nullNode?Integer.MIN_VALUE:k1.right.element.getEnd();
        k1.right.max=max(rightelement,rightleft,rightright);
        k1.right.dad=k1;
        k1.right.left.dad=k1.right;
        return k1;
    }
    
    /**
     * Rotate binary tree node with right child.
     */
    private  IntervalNode<T> rotateWithRightChild( IntervalNode<T> k1 ) {
    	//System.out.println("rotateRight");

        IntervalNode<T> k2 = k1.right;
        k2.dad=k1.dad;
        k1.right = k2.left;
        k2.left = k1;
        k2.max=k1.max;
        k2.left.max=max(k2.left.left.max,k2.left.right.max,k2.left.element.end);
        k2.left.dad=k2;
        k2.left.right.dad=k2.left;
        return k2;
    }
    private int max(int m1,int m2,int m3){
    	return m1>=m2 && m1 >=m3?m1:m2>=m3&&m2>=m1?m2:m3;
    }
    static class IntervalNode<T extends Interval<?>> {
        // Constructors
        IntervalNode( T theElement,IntervalNode<T> father ) {
        	dad=father;
            element = theElement;
            right=left=nullNode;
            level   = 1;
            if(theElement!=null)max=theElement.end;else max=Integer.MIN_VALUE;
        }
        int max;
        T element;      // The data in the node
        IntervalNode<T>     left;         // Left child
        IntervalNode<T>     right;        // Right child
        int        level;        // Level
        IntervalNode<T> dad;  //father
    }
    
     IntervalNode<T> root;
     static IntervalNode nullNode;
    static{
        nullNode = new IntervalNode( null,nullNode );
        nullNode.left = nullNode.right = nullNode;
        nullNode.level = 0;
        nullNode.max=Integer.MIN_VALUE;
    }
    private  IntervalNode<T> deletedNode;
    private  IntervalNode<T> lastNode;
    
    public void parseTree(IntervalNode<T> t){
    	if(t==nullNode)return;
    	System.out.println(t.element.start+" "+t.element.end+" "+t.max);

    	parseTree(t.left);

    	parseTree(t.right);
    }
    public static class  Flank<T>{
    	public T prev=null;
    	public T next=null;
    	public Flank(){
    	}
    	public Flank(T Prev,T Next){
    		prev=Prev;
    		next=Next;
    	}
    }
    
    // Test program; should print min and max and nothing else
    public static void main( String [ ] args ) {
        IntervalTree<Info> t = new IntervalTree<Info>( );
        t.insert(new Info(100,200,"1",'+'));
        t.insert(new Info(200,250,"2",'+'));
        t.insert(new Info(300,350,"3",'+'));
        Flank<Info> f=t.findFlank(new Info(200,250,"1",'+'));
        System.out.println(f.prev+"\n"+f.next);
        
       /* final int NUMS = 40000;
        final int GAP  =   307;
//        
//        System.out.println( "Checking... (no bad output means success)" );
//        
//        t.insert( new Interval( NUMS * 2,NUMS * 2+100 ) );
//        t.insert( new Interval( NUMS * 3,NUMS * 3+100 ) );
//        for( int i = GAP; i != 0; i = ( i + GAP ) % NUMS )
//            t.insert( new Interval( i,i+100 ) );
//        System.out.println( "Inserts complete.");
//        
//        t.remove( t.findMax( ) );
//        for( int i = 1; i < NUMS; i+= 2 )
//            t.remove( new Interval( i,i+100 ) );
//        t.remove( t.findMax( ) );
//        System.out.println( "Removes complete" );
//        
//        
//        if( ((Interval)(t.findMin( ))).start != 2 ||
//                ((Interval)(t.findMax( ))).start != NUMS - 2 )
//            System.out.println( "FindMin or FindMax error!" );
//        
//        for( int i = 2; i < NUMS; i+=2 )
//            if( ((Interval)t.find( new Interval( i,i+100 ) )).start != i )
//                System.out.println( "Error: find fails for " + i );
//        
//        for( int i = 1; i < NUMS; i+=2 )
//            if( t.find( new Interval( i,i+100 ) )  != null )
//                System.out.println( "Error: Found deleted item " + i );*/
        t.insert(new Info(1,2,"bla",'+'));
        //System.out.println("1,2");

        t.insert(new Info(2,3,"bla",'+'));
        //System.out.println("2,3");
        //t.parseTree(t.root);
        //System.out.println("_________________");
        t.insert(new Info(3,4,"bla",'+'));
        //System.out.println("3,4");
       // t.parseTree(t.root);
        //System.out.println("_________________");
        t.insert(new Info(4,5,"bla",'+'));
        //t.parseTree(t.root);
        //System.out.println("_________________");
        t.insert(new Info(5,10,"bla",'+'));
        System.out.println("_____________");

        t.insert(new Info(2,7,"bla",'+'));
        System.out.println("_____________");

        t.insert(new Info(6,8,"bla",'+'));
        System.out.println("_____________");

        t.insert(new Info(7,11,"bla",'+'));
        System.out.println("_____________");

        t.insert(new Info(2,4,"bla",'+'));
        System.out.println("_____________");

        t.insert(new Info(15,12,"bla",'+'));
        System.out.println("_____________");

        t.insert(new Info(1,9,"bla",'+'));
        System.out.println("_____________");
        t.insert(new Info(-3,8,"bla",'+'));
        t.remove(new Info(1,9,"bla",'+'));
        t.remove(new Info(-3,8,"bla",'+'));
        t.remove(new Info(15,12,"bla",'+'));
        //t.insert(new Info(-4,11,"bla"));
        //t.insert(new Info(-2,4,"bla"));
        System.out.println(t.successor(new Info(7,11,"bla",'+')));
        //t.parseTree(t.root);
    }
}








