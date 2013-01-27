import com.cycling74.max.*;


public class JittersplineExample extends MaxObject {
	
	private float px[];
	private float py[];
	private float pz[];
	private float pRed[];
	private float pGreen[];
	private float pBlue[];
	private float pScale[];
	
	private boolean lineSmooth = true;
	private boolean wireFrame = false;
	private int pointCount = 8;
	private int slices = 20;
	private int order = 3;
	private boolean outLine = false;
	private float[] outColor = {0,0,0,1};
	private boolean varyScale = true;
	private boolean varyColor = true;
	private float scaleRange = 0.2f;
	
	// morphing
	private float mx[];
	private float my[];
	private float mz[];
	private float mRed[];
	private float mGreen[];
	private float mBlue[];
	private float mScale[];
	private float slew = 0.1f;
	
	public JittersplineExample() {
		declareIO(1,1);
		px = new float[pointCount];
		py = new float[pointCount];
		pz = new float[pointCount];
		pRed = new float[pointCount];
		pBlue = new float[pointCount];
		pGreen = new float[pointCount];
		pScale = new float[pointCount];
	}
	
	public void generateRainbow() {
		post("generateRainbow()");
		for(int i=0; i<pointCount; i++) {
			px[i] = (float) Math.random()*2-1.f;
			py[i] = (float) Math.random()*2-1.f;
			pz[i] = (float) Math.random()*2-1.f;
			pRed[i] = (float) Math.random();
			pBlue[i] = (float) Math.random();
			pGreen[i] = (float) Math.random();
			pScale[i] = (float) Math.random();
		}
	}
	
	public void bang() {
		generateRainbow();
		draw();
	}

	
	public void draw() {
//		post("draw()");
		outlet(0, "reset");
		
		// "line_smooth" 
		if(lineSmooth) outlet(0,"glenable", "line_smooth");
		else outlet(0, "gldisable", "line_smooth");
		
		if(wireFrame) outlet(0,"glpolygonmode", "front_and_back line");
		else outlet(0,"glpolygonmode", "front_and_back fill");
		
		outlet(0, "beginstroke", "basic2d");
		outlet(0, "strokeparam", new Atom[]{Atom.newAtom("slices"),Atom.newAtom(slices)});
		outlet(0, "strokeparam", new Atom[]{Atom.newAtom("order"),Atom.newAtom(order)});
		outlet(0, "strokeparam", new Atom[]{Atom.newAtom("color"),Atom.newAtom(pRed[0]),Atom.newAtom(pGreen[0]),Atom.newAtom(pBlue[0])});
		outlet(0, "strokeparam", new Atom[]{Atom.newAtom("outline"),Atom.newAtom(1)});
		
		if(outLine) outlet(0, "strokeparam", new Atom[]{Atom.newAtom("outcolor"),Atom.newAtom(outColor[0]),Atom.newAtom(outColor[1]),Atom.newAtom(outColor[2]),Atom.newAtom(outColor[3])});
		outlet(0, "strokeparam", new Atom[]{Atom.newAtom("scale"),Atom.newAtom(pScale[0]*scaleRange)});
		
		for(int i=0; i<pointCount; i++) {
			if(varyScale) outlet(0, "strokeparam", new Atom[]{Atom.newAtom("scale"),Atom.newAtom(pScale[i]*scaleRange)});
			if(varyColor) outlet(0, "strokeparam", new Atom[]{Atom.newAtom("color"),Atom.newAtom(pRed[i]),Atom.newAtom(pGreen[i]),Atom.newAtom(pBlue[i])});
			outlet(0, "strokepoint", new Atom[]{Atom.newAtom(px[i]),Atom.newAtom(py[i]),Atom.newAtom(pz[i])});
		}
		
		outlet(0, "endstroke");
		outlet(0, "drawimmediate");
	}
	
	
	public void linesmooth(int v) {
	    lineSmooth = (v==1) ? true : false;
	}
	
	public void outcolor(float r, float g, float b) {
		outColor[0] = r;
		outColor[1] = g;
		outColor[2] = b;
	}
	
	public void outline(int v) {
		outLine = (v==1) ? true : false;
	}

}
