package util;

public class Shuffle {
	public static String shuffle(String seq){
		StringBuilder source=new StringBuilder(seq);
		StringBuilder target=new StringBuilder();
		
		while(source.length()>0){
			int randPos=(int)(Math.random()*source.length());
			target.append(source.charAt(randPos));
			source.deleteCharAt(randPos);
		}
		
		return target.toString();
	}
}
