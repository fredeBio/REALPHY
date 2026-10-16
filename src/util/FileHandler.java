package util;

import java.io.*;
import java.util.*;


public class FileHandler {
	
	public static void merge(File file1,File file2,File out){
		ArrayList<File> files=new ArrayList<File>();
		files.add(file1);
		files.add(file2);
		merge(files,out);
	}
	
	public static void merge(ArrayList<File> files,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			for(int i=0;i<files.size();i++){
				BufferedReader br=new BufferedReader(new FileReader(files.get(i)));
				String line="";
				while((line=br.readLine())!=null){
					bw.write(line+"\n");
				}
				br.close();
			}
			
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
    // Deletes all files and subdirectories under dir.
    // Returns true if all deletions were successful.
    // If a deletion fails, the method stops attempting to delete and returns false.
    public static boolean deleteFolder(File folder) {
        if (folder.isDirectory()) {
            String[] children = folder.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteFolder(new File(folder, children[i]));
                if (!success) {
                    return false;
                }
            }
        }
    
        // The directory is now empty so delete it
        return folder.delete();
    }
    
    public static void copy(File source, File destination){
    	try{
    		InputStream in = new FileInputStream(source);
    		OutputStream out = new FileOutputStream(destination);

    		byte[] buf = new byte[1024];
    		int len;
    		while ((len = in.read(buf)) > 0){
    			out.write(buf, 0, len);
    		}
    		in.close();
    		out.close();
    	}
    	catch(FileNotFoundException ex){
    		System.out.println(ex.getMessage() + " in the specified directory.");
    		System.exit(0);
    	}
    	catch(IOException e){
    		System.out.println(e.getMessage());  
    	}
    }

    
}
