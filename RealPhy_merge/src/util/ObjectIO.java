package util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class ObjectIO {
	public static boolean writeObject(File out,Object o){
		try{
			FileOutputStream fos=new FileOutputStream(out);
			ObjectOutputStream oos=new ObjectOutputStream(fos);
			oos.writeObject(o);
			oos.close();
			return true;
			
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}
	
	public static Object readObject(File in){
		try{
			FileInputStream fis=new FileInputStream(in);
			ObjectInputStream ois=new ObjectInputStream(fis);
			Object o= ois.readObject();
			ois.close();
			return o;
		}catch(ClassNotFoundException e){
			e.printStackTrace();
			return null;
		}catch(IOException e){
			e.printStackTrace();
			return null;
		}
	}
}
