package util;

/**
 * Exception class for failed finds/removes in search
 * trees, hash tables, and list and tree iterators.
 * @author Mark Allen Weiss
 */
public class ItemNotFoundException extends RuntimeException {
    /**
     * Construct this exception object.
     */

	private static final long serialVersionUID = 5884749257212781161L;

	public ItemNotFoundException( ) {
        super( );
    }
    
    /**
     * Construct this exception object.
     * @param message the error message.
     */
    public ItemNotFoundException( String message ) {
        super( message );
    }
}