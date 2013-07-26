package util;

import java.io.*;
import java.util.ArrayList;



import rcaller.*;

public class R_functions {
	
	public static void writeCode(RCode rc,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			bw.write(rc.getCode().toString());
			
			
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void plotTree(RCode rc,File tree){
		rc.addRCode("library(ape)");
		rc.addRCode("tree<-read.tree(\""+tree+"\")");
		rc.addRCode("plot.phylo(tree)");
	}
	
	public static void makeHistogram(RCode rc,ArrayList<Double> list,double min,double max,int bins,String name){
		ArrayList<Double> breaks=new ArrayList<Double>();
		double size=max-min;
		double steps=size/bins;
		for(int i=0;i<=bins;i++){
			breaks.add(i*steps+min);
		}
		
		rc.addDoubleArray("data", List_Array.toDouble(list));
		rc.addDoubleArray("breaks", List_Array.toDouble(breaks));
		rc.addRCode("hist(data,breaks,main=\""+name+"\")");
	}
	
	public static void main(String[] args){
		ArrayList<Double> xl=new ArrayList<Double>();
		xl.add(-5.0);
		xl.add(3.0);
		String[] test=makeTickLabelsLogDiv(makeTickPos(xl, 5),10);
		for(int i=0;i<test.length;i++){
			System.out.println(test[i]);
		}
	}
	
	public static void writeRCode(RCode rc,File out){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			String code=rc.getCode().toString();
			bw.write(code);
			bw.close();
		}catch(IOException e){
			e.printStackTrace();
		}
	}
	
    public static double[] makeTickPos(ArrayList<Double> limit,int initTickNumber){
    	int numberTicks=initTickNumber;
    	//System.out.println("limit: "+limit.get(0)+" "+limit.get(1));
    	double unrounded=Math.log10(Math.abs(limit.get(0)));
    	int startfactor=(int)(unrounded);
    	startfactor+=unrounded-startfactor<0.4?-1:0;
    	double start=limit.get(0)>0?((int)(limit.get(0)*Math.pow(10, -startfactor)))*Math.pow(10, startfactor):limit.get(0)<0?(Math.floor((limit.get(0)*Math.pow(10, -startfactor))))*Math.pow(10, startfactor):0;
    	//System.out.println(startfactor+" "+start);
    	double interval=(limit.get(1)-start)/(numberTicks);
    	long factor=Math.round(Math.log10(interval));
    	double product=interval*Math.pow(10, -factor);


    	interval=(int)product+1;

    	double step=interval*Math.pow(10,factor);

    	numberTicks=(int)((limit.get(1)-start)/step);

    	while((numberTicks-1)*step<limit.get(1)-start){
    		numberTicks++;
    	}
    	double[] list=new double[numberTicks];
    	for(int i=0;i<numberTicks;i++){
    		list[i]=(start+step*i);
    	}
    	return list;
    }
    public static String[] makeTickLabels(double[] tickPos,int roundPos){
    	String[] list=new String[tickPos.length];
    	for(int i=0;i<tickPos.length;i++){
    		list[i]=String.format("%."+roundPos+"g",tickPos[i]);
    	}
    	return list;
    }
    public static String[] makeTickLabelsEPow(double[] tickPos,int base){
    	String[] list=new String[tickPos.length];
    	for(int i=0;i<list.length;i++){
    		double logValue=tickPos[i];
    		double test=(Math.pow(base,logValue));
    		list[i]=String.format("%1$.1e",test);
    	}
    	return list;
    }
    public static String[] makeTickLabelsLogDiv(double[] tickPos,int base){
    	String[] list=new String[tickPos.length];
    	for(int i=0;i<list.length;i++){
    		double logValue=Math.abs(tickPos[i]);
    		long test=(long)(Math.pow(base,logValue));
    		list[i]=tickPos[i]<0?"1/"+test:test+"";
    	}
    	return list;
    }
    public static String[] makeTickLabelsPercent(double[] tickPos,int roundPos){
    	String[] list=new String[tickPos.length];
    	for(int i=0;i<list.length;i++){
    		double percent=tickPos[i]*100;
    		list[i]=String.format("%1$."+roundPos+"f",percent)+"%";
    	}
    	return list;
    }
    
    public static void produceScatterPlot(ArrayList<Double> l1,ArrayList<Double> l2,RCode rc,String title){
    	rc.addDoubleArray("l1", List_Array.toDouble(l1));
    	rc.addDoubleArray("l2", List_Array.toDouble(l2));
    	rc.addRCode("plot(l1,l2,main="+title+")");
    }
    
	public static void addMyImagePlot(RCode rcode,int scale){
		
		rcode.addRCode("myImagePlot <- function(x, ...){");
		rcode.addRCode(" min <- min(x)");
		rcode.addRCode("max <- max(x)");
		rcode.addRCode("if(min==max){min<-min-0.1}");
		rcode.addRCode(" yLabels <- rownames(x)");
		rcode.addRCode("  xLabels <- colnames(x)");
		rcode.addRCode(" title <-c()");
		rcode.addRCode(" layout<-TRUE");
		rcode.addRCode("  colorscale<-TRUE");
		rcode.addRCode(" # check for additional function arguments");
		rcode.addRCode(" if( length(list(...)) ){");
		rcode.addRCode("  Lst <- list(...)");
		rcode.addRCode("  if( !is.null(Lst$zlim) ){");
		rcode.addRCode("   min <- Lst$zlim[1]");
		rcode.addRCode("   max <- Lst$zlim[2]");
		rcode.addRCode("}");
		rcode.addRCode(" if( !is.null(Lst$yLabels) ){");
		rcode.addRCode("   yLabels <- c(Lst$yLabels)");
		rcode.addRCode("}");
		rcode.addRCode(" if( !is.null(Lst$xLabels) ){");
		rcode.addRCode("    xLabels <- c(Lst$xLabels)");
		rcode.addRCode(" }");
		rcode.addRCode(" if( !is.null(Lst$title) ){");
		rcode.addRCode("   title <- Lst$title");
		rcode.addRCode(" }");
		rcode.addRCode("if(!is.null(Lst$layout)){");
		rcode.addRCode("   layout<-Lst$layout");
		rcode.addRCode("}");
		rcode.addRCode(" if(!is.null(Lst$colorscale)){");
		rcode.addRCode("    colorscale<-Lst$colorscale");
		rcode.addRCode(" }");
		rcode.addRCode("}");
		rcode.addRCode("# check for null values");
		rcode.addRCode("if( is.null(xLabels) ){");
		rcode.addRCode("   xLabels <- c(1:ncol(x))");
		rcode.addRCode("}");
		rcode.addRCode("if( is.null(yLabels) ){");
		rcode.addRCode("yLabels <- c(1:nrow(x))");
		rcode.addRCode("}");

		rcode.addRCode("if(layout)layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))");
		rcode.addRCode("# Red and green range from 0 to 1 while Blue ranges from 1 to 0");
		if(scale==0){
			rcode.addRCode("ColorRamp <- rgb( c(seq(1,1,length=256)),  # Red");
			rcode.addRCode("c(seq(1,0,length=256)),  # Green");
			rcode.addRCode("c(   seq(1,0,length=256)))  # Blue");
		}else if(scale==1){
			rcode.addRCode("ColorRamp <- rgb( c(seq(1,1,length=256),seq(0,1,length=256)),  # Red");
			rcode.addRCode("c(seq(1,0,length=256),seq(1,1,length=256)),  # Green");
			rcode.addRCode("c(   seq(1,0,length=256),seq(0,1,length=256)))  # Blue");
		}else if(scale==2){
			rcode.addRCode("ColorRamp <- rgb( c(seq(1,1,length=256),seq(1,0,length=256)),  # Red");
			rcode.addRCode("c(seq(0,1,length=256),seq(1,1,length=256)),  # Green");
			rcode.addRCode("c(   seq(0,1,length=256),seq(1,0,length=256)))  # Blue");
		}else if(scale==3){
			rcode.addRCode("ColorRamp <- rgb( c(seq(1,1,length=256)),  # Red");
			rcode.addRCode("c(seq(0,1,length=256)),  # Green");
			rcode.addRCode("c(   seq(0,1,length=256)))  # Blue");
		}else if(scale==4){
			rcode.addRCode("ColorRamp <- rgb( c(seq(1,1,length=256),seq(1,0,length=256)),  # Red");
			rcode.addRCode("c(seq(0,0,length=256),seq(1,1,length=256)),  # Green");
			rcode.addRCode("c(   seq(0,0,length=256),seq(1,0,length=256)))  # Blue");
		}
		rcode.addRCode("ColorLevels <- seq(min, max, length=length(ColorRamp))");
		

		rcode.addRCode("# Data Map");
		rcode.addRCode("par(mar = c(3,5,3.5,2),cex=0.3)");
		rcode.addRCode("image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab=\"\",");
		rcode.addRCode(" ylab=\"\", axes=FALSE, zlim=c(min,max))");
		rcode.addRCode("#print text into cells");
		rcode.addRCode("posx=c()");
		rcode.addRCode("posy=c()");
		rcode.addRCode("label=c()");

		rcode.addRCode("for(i in 1:dim(x)[1]){");
		rcode.addRCode("  for(j in 1:dim(x)[2]){");
		rcode.addRCode(" posx=c(posx,i)");
		rcode.addRCode(" posy=c(posy,j)");
		rcode.addRCode(" label=c(label,x[j,i])");

		rcode.addRCode(" }");
		rcode.addRCode(" }");
		rcode.addRCode("	 text(posx,posy,labels=label)");

		rcode.addRCode("if( !is.null(title) ){");
		rcode.addRCode("   title(main=title)");
		rcode.addRCode(" }");
		rcode.addRCode(" axis(BELOW<-1, at=1:length(xLabels), labels=xLabels)");
		rcode.addRCode(" axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1)");

		rcode.addRCode(" # Color Scale");
		rcode.addRCode(" if(colorscale==TRUE){");
		rcode.addRCode("       drawColorScale(min,max)");
		rcode.addRCode(" }");
		rcode.addRCode("if(layout) layout(1)");
		rcode.addRCode("	}");
	}
	public static void runRCode(RCode rc,File scriptPath){
    	rc.addRCode("dev.off()");
    	RCaller rca=new RCaller();
    	rca.setRCode(rc);
    	rca.setRscriptExecutable(scriptPath.toString());
    	rca.runOnly();
    }
	  public static RCode plot_InitJpg(File out,int row,int col){//5,5.5,0.5
	    	RCode rc=new RCode();
	    	rc.addRCode("jpeg(file=\""+out+"\",width="+col*8+",height="+8*row+",units='cm',res=300)");
	    	rc.addRCode("par(mfcol=c("+row+","+col+"),mar=c(5, 7, 4, 2)+0.1,cex.lab=1.4,cex.axis=1.4)");
	    	rc.addRCode("library(gplots)");
	    	return rc;
	    }
	  public static RCode plot_InitJpg(File out,double height,double width){//5,5.5,0.5
	    	RCode rc=new RCode();
	    	rc.addRCode("jpeg(file=\""+out+"\",width="+width+",height="+height+",units='cm',res=300)");
	    	rc.addRCode("library(gplots)");
	    	return rc;
	    }
	  public static RCode plot_InitJpg(File out){//5,5.5,0.5
	    	RCode rc=new RCode();
	    	rc.addRCode("jpeg(file=\""+out+"\")");
	    	rc.addRCode("library(gplots)");
	    	return rc;
	    }
	  public static RCode plot_InitTiff(File out,double height,double width){//5,5.5,0.5
	    	RCode rc=new RCode();
	    	rc.addRCode("tiff(file=\""+out+"\",width="+width+",height="+height+",units='cm',res=300)");
	    	rc.addRCode("library(gplots)");
	    	return rc;
	    }
    public static void produceRcode_Boxplot(RCode rc,ArrayList<ArrayList<Double>> y,ArrayList<Double> at,ArrayList<Double> xlim,ArrayList<Double> ylim,String ylabel){
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	for(int i=0;i<y.size();i++){
    		String dsName="y"+i;
    		rc.addDoubleArray(dsName, List_Array.toDouble(y.get(i)));
    		sb.append(dsName+",");
    	}
    	rc.addDoubleArray("at", List_Array.toDouble(at));
    	rc.addDoubleArray("xl", List_Array.toDouble(xlim));
    	rc.addDoubleArray("yl", List_Array.toDouble(ylim));


    	rc.addRCode(sb+"axes=FALSE,at=at,ylab=\""+ylabel+"\", xlab=\"Divergence\",xlim=xl,ylim=yl,boxwex=0.02)");
    	rc.addRCode("axis(1)");
    	rc.addRCode("axis(2)");
    }

    public static void produceRcode_BoxplotCategories(RCode rc,ArrayList<ArrayList<Double>> y,ArrayList<String> cats,String ylabel,double ymin,double ymax){
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	for(int i=0;i<y.size();i++){
    		String dsName="y"+i;
    		rc.addDoubleArray(dsName, List_Array.toDouble(y.get(i)));
    		sb.append(dsName+",");
    	}
    	rc.addStringArray("names", List_Array.toString(cats));


    	sb.append("axes=FALSE,ylab=\""+ylabel+"\",range=0,ylim=c("+ymin+","+ymax+")");
    	rc.addRCode(sb.toString()+")");
    	//rc.addRCode("rect(0,0,"+16+","+ymax+",col=\"red\",border=NA)");
    	//rc.addRCode("rect(0,"+ymin+",16,0,col=\"green\",border=NA)");
    	//rc.addRCode(sb.toString()+",add=TRUE)");
    	rc.addRCode("abline(0,0)");
    	rc.addRCode("par(las=3)");
    	rc.addRCode("axis(1,labels=names,at=1:15)");
    	rc.addRCode("par(las=1)");
    	rc.addRCode("axis(2)");
    }
    
    private static void addDoubleArrays(ArrayList<ArrayList<Double>> list,RCode rc,StringBuffer sb,String name){
    	for(int i=0;i<list.size();i++){
    		String dsName=name+i;
    		rc.addDoubleArray(dsName, List_Array.toDouble(list.get(i)));
    		sb.append(dsName+",");
    	}
    }
    
    public static void produceRcode_BoxplotCategories(RCode rc,ArrayList<ArrayList<Double>> list1,ArrayList<ArrayList<Double>> list2,ArrayList<String> cats,String ylabel,String col1,String col2){
    	produceRcode_BoxplotCategories(rc, list1, list2, cats, ylabel, col1, col2,"");
    }
    
    public static void produceRcode_BoxplotCategories(RCode rc,ArrayList<ArrayList<Double>> list1,ArrayList<ArrayList<Double>> list2,ArrayList<String> cats,String ylabel,String col1,String col2,String ylim){
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(list1, rc, sb, "list1");
    	
    	rc.addStringArray("names", List_Array.toString(cats));

    	
    	rc.addRCode(sb+"axes=FALSE,ylab=\""+ylabel+"\",range=0,border=\""+col1+"\""+ylim+")");
    	sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(list2, rc, sb, "list2");
    	rc.addRCode(sb+"axes=FALSE,ylab=\""+ylabel+"\",range=0,add=TRUE,border=\""+col2+"\""+ylim+")");
    	rc.addRCode("abline(0,0)");
    	rc.addRCode("par(las=3)");
    	rc.addRCode("axis(1,labels=names,at=1:15)");
    	rc.addRCode("par(las=1)");
    	rc.addRCode("axis(2)");
    }
    
    public static void produceRcode_BoxplotCategories(RCode rc,ArrayList<ArrayList<Double>> list1,ArrayList<ArrayList<Double>> list2,ArrayList<String> cats,String ylabel,String col1,String col2,double ymin,double ymax){
    	produceRcode_BoxplotCategories(rc, list1, list2, cats, ylabel, col1, col2, ",ylim=c("+ymin+","+ymax+")");
    }

    private static ArrayList<ArrayList<Double>> getList(ArrayList<ArrayList<LogLikeDiff>> list,int var){
    	ArrayList<ArrayList<Double>> newList=new ArrayList<ArrayList<Double>>();
    	for(int i=0;i<list.size();i++){
    		ArrayList<Double> temp=new ArrayList<Double>();
    		for(int j=0;j<list.get(i).size();j++){
    			temp.add(list.get(i).get(j).get(var));
    		}
    		newList.add(temp);
    	}
    	return newList;
    		
    }
    
    private static ArrayList<ArrayList<Double>> getList(ArrayList<ArrayList<LogLikeDiff>> list,int var,ArrayList<Integer> items){
    	ArrayList<ArrayList<Double>> newList=new ArrayList<ArrayList<Double>>();
    	for(int i=0;i<list.size();i++){
    		ArrayList<Double> temp=new ArrayList<Double>();
    		for(int j=0;j<items.size();j++){
    			temp.add(list.get(i).get(items.get(j)).get(var));
    		}
    		newList.add(temp);
    	}
    	return newList;
    		
    }
    
    public static void produceRcode_BoxplotCategoriesLogLike(RCode rc,ArrayList<ArrayList<LogLikeDiff>> list,ArrayList<String> cats,String ylabel,String col1,String col2,double yminDiff,double ymaxDiff,double yminRP,double ymaxRP,double yminExact,double ymaxExact){
    	ArrayList<ArrayList<Double>> fullCorr1=getList(list,LogLikeDiff.FC);
    	ArrayList<ArrayList<Double>> fullInc1=getList(list,LogLikeDiff.FI);
    	ArrayList<ArrayList<Double>> rpCorr1=getList(list,LogLikeDiff.RC);
    	ArrayList<ArrayList<Double>> rpInc1=getList(list,LogLikeDiff.RI);
    	ArrayList<ArrayList<Double>> total1=getList(list,LogLikeDiff.TO);
    	
    	rc.addRCode("par(mfcol=c(3,1))");
    	produceRcode_BoxplotCategories(rc, fullCorr1, fullInc1, cats, ylabel, col1, col2,yminExact,ymaxExact);
    	produceRcode_BoxplotCategories(rc, rpCorr1, rpInc1, cats, ylabel, col1, col2,yminRP,ymaxRP);
    	produceRcode_BoxplotCategories(rc, total1, cats, ylabel,yminDiff,ymaxDiff);
    	
    }
    
    public static void produceRcode_BoxplotCategoriesLogLike(RCode rc,ArrayList<ArrayList<LogLikeDiff>> list,ArrayList<Integer> items1,ArrayList<Integer> items2,ArrayList<String> cats,String ylabel,String col1,String col2,double yminDiff,double ymaxDiff,double yminRP,double ymaxRP,double yminExact,double ymaxExact){
    	ArrayList<ArrayList<Double>> fullCorr1=getList(list,LogLikeDiff.FC,items1);
    	
    	ArrayList<ArrayList<Double>> fullInc1=getList(list,LogLikeDiff.FI,items1);
    	ArrayList<ArrayList<Double>> rpCorr1=getList(list,LogLikeDiff.RC,items1);
    	ArrayList<ArrayList<Double>> rpInc1=getList(list,LogLikeDiff.RI,items1);
    	ArrayList<ArrayList<Double>> total1=getList(list,LogLikeDiff.TO,items1);

    	ArrayList<ArrayList<Double>> fullCorr2=getList(list,LogLikeDiff.FC,items2);
    	ArrayList<ArrayList<Double>> fullInc2=getList(list,LogLikeDiff.FI,items2);
    	ArrayList<ArrayList<Double>> rpCorr2=getList(list,LogLikeDiff.RC,items2);
    	ArrayList<ArrayList<Double>> rpInc2=getList(list,LogLikeDiff.RI,items2);
    	ArrayList<ArrayList<Double>> total2=getList(list,LogLikeDiff.TO,items2);

    	
    	rc.addRCode("par(mfcol=c(3,2))");
    	produceRcode_BoxplotCategories(rc, fullCorr1, fullInc1, cats, ylabel, col1, col2,yminExact,ymaxExact);
    	produceRcode_BoxplotCategories(rc, rpCorr1, rpInc1, cats, ylabel, col1, col2,yminRP,ymaxRP);
    	produceRcode_BoxplotCategories(rc, total1, cats, ylabel,yminDiff,ymaxDiff);
    	
    	produceRcode_BoxplotCategories(rc, fullCorr2, fullInc2, cats, ylabel, col1, col2,yminExact,ymaxExact);
    	produceRcode_BoxplotCategories(rc, rpCorr2, rpInc2, cats, ylabel, col1, col2,yminRP,ymaxRP);
    	produceRcode_BoxplotCategories(rc, total2, cats, ylabel,yminDiff,ymaxDiff);

    }
    
    public static void produceRcode_MLNI(RCode rc,ArrayList<Double> x,ArrayList<Double> y,ArrayList<Double> xl,ArrayList<Double> yl,String ylabel,char pch){
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("y", List_Array.toDouble(y));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	rc.addRCode("plot(x,y,pch=\""+pch+"\",type=\"p\",col=1,ylim=yl,xlim=xl,axes=FALSE,xlab=\"Divergence\",ylab=\""+ylabel+"\")");
    	rc.addRCode("axis(1)");
    	rc.addRCode("axis(2)");
    }
    
	
    public static void produceRcode_Binned(RCode rc,ArrayList<Double> x,ArrayList<Double> y,ArrayList<Double> stderr,ArrayList<Double> xl,ArrayList<Double> yl,String ylabel){
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("y", List_Array.toDouble(y));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	rc.addDoubleArray("stderr", List_Array.toDouble(stderr));
    	rc.addRCode("plotCI(x,y,uiw=stderr,gap=0,err=\"y\",ylim=yl,ylab=\""+ylabel+"\", xlab=\"Divergence\",axes=FALSE,xlim=xl)");
    	rc.addRCode("axis(1)");
    	rc.addRCode("axis(2)");
    }
	
	public static void produceRcode_Binned(RCode rc,ArrayList<Double> x,ArrayList<Double> y,ArrayList<Double> stderr,ArrayList<Double> xMP,ArrayList<Double> yMP,ArrayList<Double> stderrMP,ArrayList<Double> xl,ArrayList<Double> yl,String ylabel){
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("y", List_Array.toDouble(y));
    	rc.addDoubleArray("xMP", List_Array.toDouble(xMP));
    	rc.addDoubleArray("yMP", List_Array.toDouble(yMP));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	rc.addDoubleArray("stderr", List_Array.toDouble(stderr));
    	rc.addDoubleArray("stderrMP", List_Array.toDouble(stderrMP));
    	rc.addRCode("plotCI(x,y,uiw=stderr,gap=0,err=\"y\",ylim=yl,ylab=\""+ylabel+"\", xlab=\"Divergence\",axes=FALSE,xlim=xl)");
    	rc.addRCode("plotCI(xMP,yMP,gap=0,uiw=stderrMP,add=TRUE,col=\"red\")");
    	rc.addRCode("axis(1)");
    	rc.addRCode("axis(2)");
    }

    public static void produceRcode_BoxplotBinned(RCode rc,ArrayList<Double> x,ArrayList<ArrayList<Double>> y,ArrayList<Double> xl,ArrayList<Double> yl,String ylabel,String title,String col1){
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(y, rc, sb, "list1");
    	
    	rc.addDoubleArray("x", List_Array.toDouble(x));

    	
    	rc.addRCode(sb+"axes=FALSE,ylab=\""+ylabel+"\",range=0,border=\""+col1+"\",ylim=yl,xlim=xl,at=x,boxwex=0.02,main=\""+title+"\")");
    	rc.addRCode("par(las=1)");
    	rc.addRCode("axis(1)");
    	rc.addRCode("par(las=1)");
    	rc.addRCode("axis(2)");
    }
    public static void produceRcode_BoxplotBinned(RCode rc,ArrayList<Double> x,ArrayList<ArrayList<Double>> y,ArrayList<Double> xl,ArrayList<Double> yl,String title,String xAxis,String yAxis,String xTickLabels[],String yTickLabels[],double[] xTickPos,double[] yTickPos,String colour,boolean axesLabels,double xlabPos,double ylabPos){
    	if(title.length()==0)rc.addRCode("par(mar=c(4.5, 9, 1, 1)+0.1,cex.lab=2,cex.axis=2)");
    	else rc.addRCode("par(mar=c(4.5, 9, 2, 1)+0.1,cex.lab=2,cex.axis=2,cex.main=2)");
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	rc.addDoubleArray("xTickPos",xTickPos);
    	rc.addDoubleArray("yTickPos", yTickPos);
    	rc.addStringArray("xTickLabels", xTickLabels);
    	rc.addStringArray("yTickLabels",yTickLabels);
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(y, rc, sb, "list1");
    	
    	rc.addDoubleArray("x", List_Array.toDouble(x));

    	rc.addRCode(sb+"axes=FALSE,range=0,ylim=yl,xlim=xl,at=x,main=\""+title+"\",boxwex=0.02,border=\""+colour+"\")");
    	if(axesLabels){
    		
    		rc.addRCode("par(las=1)");
    		rc.addRCode("title(ylab=\""+yAxis+"\",mgp=c("+ylabPos+",1,0))");
    		rc.addRCode("axis(1,at=xTickPos,labels=xTickLabels)");
    		rc.addRCode("title(xlab=\""+xAxis+"\",mgp=c("+xlabPos+",1,0))");
    		rc.addRCode("par(las=1)");
    		rc.addRCode("axis(2,at=yTickPos,labels=yTickLabels)");
    	}else{
    		rc.addRCode("axis(1,at=xTickPos,labels=FALSE)");
    		rc.addRCode("axis(2,at=yTickPos,labels=FALSE)");

    	}
    }
    public static void produceRcode_BoxplotBinned(RCode rc,ArrayList<Double> x,ArrayList<ArrayList<Double>> y,ArrayList<Double> xMP,ArrayList<ArrayList<Double>> yMP,ArrayList<Double> xl,ArrayList<Double> yl,String col1,String col2,String xAxis,String yAxis,String xTickLabels[],String yTickLabels[],double[] xTickPos,double[] yTickPos,boolean axesLabels){
    	rc.addRCode("par(cex.axis=2,cex.lab=2,mar=c(5,7,1,1)+0.1)");
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("xMP", List_Array.toDouble(xMP));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	rc.addDoubleArray("xTickPos",xTickPos);
    	rc.addDoubleArray("yTickPos", yTickPos);
    	rc.addStringArray("xTickLabels", xTickLabels);
    	rc.addStringArray("yTickLabels",yTickLabels);
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(y, rc, sb, "list1");
    	
    	rc.addDoubleArray("x", List_Array.toDouble(x));

    	
    	rc.addRCode(sb+"axes=FALSE,range=0,border=\""+col1+"\",ylim=yl,xlim=xl,at=x,boxwex=0.02)");
    	sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(yMP, rc, sb, "list2");
    	rc.addRCode(sb+"axes=FALSE,range=0,add=TRUE,border=\""+col2+"\",ylim=yl,xlim=xl,at=xMP,boxwex=0.02)");
    	if(axesLabels){
    		rc.addRCode("par(las=1)");
    		rc.addRCode("title(xlab=\""+xAxis+"\",mgp=c(3,1,0))");
    		rc.addRCode("title(ylab=\""+yAxis+"\",mgp=c(4.5,1,0))");

    		rc.addRCode("axis(1,at=xTickPos,labels=xTickLabels)");
    		rc.addRCode("par(las=1)");
    		rc.addRCode("axis(2,at=yTickPos,labels=yTickLabels)");
    	}else{
    		rc.addRCode("axis(1,at=xTickPos)");
    		rc.addRCode("axis(2,at=yTickPos)");

    	}
    }
    public static void produceRcode_BoxplotBinned(RCode rc,ArrayList<Double> x,ArrayList<ArrayList<Double>> y,ArrayList<Double> xMP,ArrayList<ArrayList<Double>> yMP,ArrayList<Double> xl,ArrayList<Double> yl,String ylabel,String col1,String col2){
    	rc.addDoubleArray("x", List_Array.toDouble(x));
    	rc.addDoubleArray("xMP", List_Array.toDouble(xMP));
    	rc.addDoubleArray("xl", List_Array.toDouble(xl));
    	rc.addDoubleArray("yl", List_Array.toDouble(yl));
    	StringBuffer sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(y, rc, sb, "list1");
    	
    	rc.addDoubleArray("x", List_Array.toDouble(x));

    	
    	rc.addRCode(sb+"axes=FALSE,ylab=\""+ylabel+"\",range=0,border=\""+col1+"\",ylim=yl,xlim=xl,at=x,boxwex=0.02)");
    	sb=new StringBuffer();
    	sb.append("boxplot(");
    	addDoubleArrays(yMP, rc, sb, "list2");
    	rc.addRCode(sb+"axes=FALSE,ylab=\""+ylabel+"\",range=0,add=TRUE,border=\""+col2+"\",ylim=yl,xlim=xl,at=xMP,boxwex=0.02)");
    	//rc.addRCode("abline(0,0)");
    	rc.addRCode("par(las=1)");
    	rc.addRCode("axis(1)");
    	rc.addRCode("par(las=1)");
    	rc.addRCode("axis(2)");
    }
	
}
