package util;

public class WordFrequency implements Comparable<WordFrequency> {
	public int frequency;
	public String word;
	
	public WordFrequency(int freq,String Word){
		word=Word;
		frequency=freq;
	}
	
	public int compareTo(WordFrequency o) {
		// TODO Auto-generated method stub
		if(o.frequency>frequency){
			return 1;
		}else if(o.frequency<frequency){
			return -1;
		}else if(o.frequency==frequency){
			return word.compareTo(o.word);
		}
		return 0;
	}
	@Override
	public String toString(){
		return word+"\t"+frequency;
	}

}
