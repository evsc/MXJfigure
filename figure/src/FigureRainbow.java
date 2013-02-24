/*
 * Figure de la terre - mxj objects
 * 
 * make rainbow curves 
 * with jitter jit.gl.sketch stroke function
 * 
 * simple version: no animation, but morphing
 * TODO update math and drawing similar to FigureRainbowBezierTriple
 * 
 */


import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;


public class FigureRainbow extends MaxObject {
	
	// render context for opengl objects
	private String context = "foo";
	
	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	
	
	// position, color and scale arrays for the main bezier curve
	private float px[];
	private float py[];
	private float pz[];
	private float pRed[];
	private float pGreen[];
	private float pBlue[];
	private float pScale[];
	
	// display parameters
	private boolean lineSmooth = true;
	private boolean wireFrame = false;
	private int pointCount = 3;
	private int maxPoints = 10;
	private int bSlices = 20;		// bezier slices
	private int bOrder = 3;			// bezier order (1st order, 2nd order, ... bezier curves)
	private boolean outLine = false;
	private float[] outColor = {0,0,0,1};
	
	
	private boolean varyScale = false;
	private boolean varyColor = false;
	private float scaleRange = 0.1f;
	private boolean showPoints = false;
	private boolean colorMode = true;
	
	// rainbow
	private int rBands = 7;
	private float rWidth = 0.2f;
	private float[] rBase = {1,0,0};	// direction of bands distribution on base
	private float[] rTop = {0,-1,0};		// direction of bands distribution on arch
	private float[] rPosition = {0,0,0};	// center position of rainbow (arch center on ground)
	private float[] rSize = {1,1,1};		// scale of rainbow
	private float rDir = 1;					// direction. -1 inverts direction of bands
	private int rMode = 1;			// rainbow building mode
	
	// position and scale arrays for morphing
	private boolean morphing = false;
	private float slew = 0.1f;
	private float mx[];
	private float my[];
	private float mz[];
	private float mScale[];
	private float mPosition[] = {0,0,0};
	private float mSize[] = {1,1,1};
	private float mWidth = 0.2f;
	private float mBands = 7;
	private float mScaleRange = 0.1f;
	
	/* instantiating mxj without argument causes a problem */
	public FigureRainbow() {
		// bail method presents instantiation of object and prints error message 
		bail("please provide render context as argument");
	}
	
	/* instantiate mxj with render context as argument */
	public FigureRainbow(String c) {
		context = c;			// 	render context
		declareIO(2,1);			/* 	delare number of inlets and outlets, 
 									of DataTypes.ALL */
		
		// assist message for inlets and outlets (for mouse hover)
		setInletAssist(new String[] {"bang to compute and draw", "input settings"});
		setOutletAssist(new String[] {"outputs jit.gl.texture object"});
		
		
		// instantiate jitter objects
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", 1);
		sketch.setAttr("glclearcolor", new Atom[]{Atom.newAtom(0.),Atom.newAtom(0.),Atom.newAtom(0.),Atom.newAtom(1.)});
		sketch.setAttr("fsaa", 1);
		sketch.send("automatic", 0);	/* set to not-automatic, to be able to use
										   begin_capture and drawimmediate */

		texture = new JitterObject("jit.gl.texture");
		texture.setAttr("drawto", context);
//		texture.setAttr("name", "tex_rainbow");	// don't set name, to avoid doubles
		texture.setAttr("dim", new Atom[]{Atom.newAtom(640),Atom.newAtom(480)});
		
		// set arrays to max number of pointcount
		int m = maxPoints;
		px = new float[m];
		py = new float[m];
		pz = new float[m];
		pRed = new float[m];
		pBlue = new float[m];
		pGreen = new float[m];
		pScale = new float[m];
		mx = new float[m];
		my = new float[m];
		mz = new float[m];
		mScale = new float[m];
	}
	
	/* defines the points of the rainbow randomly */
	public void generateRainbow() {
		post("generateRainbow() mode:"+rMode);
		
		switch(rMode) {
		case 0:	// random points
			for(int i=0; i<pointCount; i++) {
				px[i] = (float) Math.random()*2-1.f;
				py[i] = (float) Math.random()*2-1.f;
				pz[i] = (float) Math.random()*2-1.f;
				pScale[i] = (float) Math.random();
			}
			break;
		case 1: // strict rainbow from 4 points
				pointCount = 4;
				px[0] = px[1] = rPosition[0] - rSize[0];
				px[2] = px[3] = rPosition[0] + rSize[0];
				py[0] = py[3] = rPosition[1];
				py[1] = py[2] = rPosition[1] + rSize[1];
				pz[0] = pz[1] = pz[2] = pz[3] = 0;
				for(int i=0; i<pointCount; i++) pScale[i] = 1.0f;
			break;
			
		case 2: // strict rainbow from 3 points
				pointCount = 3;
				px[0] = rPosition[0] - rSize[0];
				px[2] = rPosition[0] + rSize[0];
				px[1] = rPosition[0];
				py[0] = py[2] = rPosition[1];
				py[1] = rPosition[1] + rSize[1];
				pz[0] = pz[1] = pz[2] = 0;
				for(int i=0; i<pointCount; i++) pScale[i] = 1.0f;
			break;
		
		case 3: // advance array of points, one new random point
				for(int i=0; i<pointCount; i++) {
					if(i<pointCount-1) {
						px[i] = px[i+1];
						py[i] = py[i+1];
						pz[i] = pz[i+1];
					}
					pScale[i] = 1.0f;
				}
				px[pointCount-1] = (float) Math.random()*2-1.f;
				py[pointCount-1] = (float) Math.random()*2-1.f;
				pz[pointCount-1] = (float) Math.random()*2-1.f;
			break;
		
			
		}
	}

	
	/* draws and captures the rainbow drawing to sketch object, and outputs sketch as jitter texture object */
	public void bang() {
		texture.call("begin_capture");	// begin capturing	
		draw();							// draw rainbow to sketch
		texture.call("end_capture");			// end capturing
		texture.call("draw");					// to output texture? 
		outlet(0,"jit_gl_texture",texture.getAttr("name"));		// output texture
	}

	/* draw rainbow to jitter sketch object */
	public void draw() {
		
		// create local variables, to avoid conflict when life-updating variables while rendering
		float _px[] = new float[maxPoints];
		float _py[] = new float[maxPoints];
		float _pz[] = new float[maxPoints];
		
		for(int i=0; i<pointCount; i++) {
			mx[i] = (morphing) ? mx[i] += (px[i]-mx[i])*slew : px[i];
			my[i] = (morphing) ? my[i] += (py[i]-my[i])*slew : py[i];
			mz[i] = (morphing) ? mz[i] += (pz[i]-mz[i])*slew : pz[i];
			_px[i] = mx[i];
			_py[i] = my[i];
			_pz[i] = mz[i];
		}
		
		mBands = (morphing) ? mBands += (rBands - mBands)*slew : rBands;
		int _bands = (int) Math.floor(mBands);
		mWidth = (morphing) ? mWidth += (rWidth - mWidth)*slew : rWidth;
		float _w = mWidth;
		mScaleRange = (morphing) ? mScaleRange += (scaleRange - mScaleRange)*slew : scaleRange;
		float _sr = mScaleRange;
		
		mPosition[0] = (morphing) ? mPosition[0] += (rPosition[0] - mPosition[0])*slew : rPosition[0];
		mPosition[1] = (morphing) ? mPosition[1] += (rPosition[1] - mPosition[1])*slew : rPosition[1];
		mPosition[2] = (morphing) ? mPosition[2] += (rPosition[2] - mPosition[2])*slew : rPosition[2];
		float _x = mPosition[0];
		float _y = mPosition[1];
		float _z = mPosition[2];
		mSize[0] = (morphing) ? mSize[0] += (rSize[0] - mSize[0])*slew : rSize[0];
		mSize[1] = (morphing) ? mSize[1] += (rSize[1] - mSize[1])*slew : rSize[1];
		mSize[2] = (morphing) ? mSize[2] += (rSize[2] - mSize[2])*slew : rSize[2];
		float _sx = mSize[0];
		float _sy = mSize[1];
		float _sz = mSize[2];
		
		float _d = rDir;
		
		
		sketch.call("reset");
		
		// "line_smooth" 
		if(lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");
		
		if(wireFrame) sketch.call("glpolygonmode", new Atom[]{Atom.newAtom("front_and_back"),Atom.newAtom("line")});
		else sketch.call("glpolygonmode", new Atom[]{Atom.newAtom("front_and_back"),Atom.newAtom("fill")});
		
		
		
		for(int b=0; b<_bands; b++) {
			
			float index = b * 1.0f / (float) (_bands-1);
			float cindex =  b * 1.0f / (float) _bands;
			
			// color
			Atom[] bandColor = new Atom[]{Atom.newAtom("color"),Atom.newAtom(rRed(cindex)),Atom.newAtom(rGreen(cindex)),Atom.newAtom(rBlue(cindex))};
			
		
			sketch.call("beginstroke", "basic2d");
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("slices"),Atom.newAtom(bSlices)});
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("order"),Atom.newAtom(bOrder)});
			sketch.call("strokeparam", bandColor);
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("outline"),Atom.newAtom(1)});
			
			if(outLine) sketch.call("strokeparam", new Atom[]{Atom.newAtom("outcolor"),Atom.newAtom(outColor[0]),Atom.newAtom(outColor[1]),Atom.newAtom(outColor[2]),Atom.newAtom(outColor[3])});
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("scale"),Atom.newAtom(pScale[0]*_sr)});
			
			// start of the rainbow. first point
			if(varyScale) sketch.call("strokeparam", new Atom[]{Atom.newAtom("scale"),Atom.newAtom(pScale[0]*_sr)});
			sketch.call("strokepoint", new Atom[]{Atom.newAtom(_x+_px[0]*_sx+rBase[0]*index*_w*_d),Atom.newAtom(_y+_py[0]*_sy+rBase[1]*index*_w*_d),Atom.newAtom(_z+_pz[0]*_sz+rBase[2]*index*_w*_d)});
			
			// inbetween of the rainbow
			for(int i=1; i<pointCount-1; i++) {
				if(varyScale) sketch.call("strokeparam", new Atom[]{Atom.newAtom("scale"),Atom.newAtom(pScale[i]*_sr)});
				float _r = (_px[i] > _x) ? -1.f : 1.f;	// add x-base on points on left, subtract on points on right
//				if(varyColor) sketch.call("strokeparam", new Atom[]{Atom.newAtom("color"),Atom.newAtom(pRed[i]),Atom.newAtom(pGreen[i]),Atom.newAtom(pBlue[i])});
				sketch.call("strokepoint", new Atom[]{Atom.newAtom(_x+_px[i]*_sx+rTop[0]*index*_w*_d + _r*rBase[0]*index*_w*_d),Atom.newAtom(_y+_py[i]*_sy+rTop[1]*index*_w*_d),Atom.newAtom(_z+_pz[i]*_sz+rTop[2]*index*_w*_d)});
			}
			
			int i = pointCount-1;
			// end of the rainbow. last point. subtract rBase value
			if(varyScale) sketch.call("strokeparam", new Atom[]{Atom.newAtom("scale"),Atom.newAtom(pScale[i]*_sr)});
			sketch.call("strokepoint", new Atom[]{Atom.newAtom(_x+_px[i]*_sx-rBase[0]*index*_w*_d),Atom.newAtom(_y+_py[i]*_sy-rBase[1]*index*_w*_d),Atom.newAtom(_z+_pz[i]*_sz-rBase[2]*index*_w*_d)});
			
			sketch.call("endstroke");
		
		}
		
		// display the control points of the main bezier curve
		if(showPoints) {
			sketch.call("beginstroke", "line");
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("order"),Atom.newAtom(1)});
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("segments"),Atom.newAtom(1)});
			sketch.call("strokeparam", new Atom[]{Atom.newAtom("color"),Atom.newAtom(1),Atom.newAtom(1),Atom.newAtom(1)});
			for(int i=0; i<pointCount; i++) {
				sketch.call("strokepoint", new Atom[]{Atom.newAtom(_x+_px[i]*_sx),Atom.newAtom(_y+_py[i]*_sy),Atom.newAtom(_z+_pz[i]*_sz)});
			}
			sketch.call("endstroke");
			for(int i=0; i<pointCount; i++) {
				sketch.call("moveto", new Atom[]{Atom.newAtom(_x+_px[i]*_sx),Atom.newAtom(_y+_py[i]*_sy),Atom.newAtom(_z+_pz[i]*_sz)});
				sketch.call("circle", new Atom[]{Atom.newAtom(0.03)});
			}
		}
		
		sketch.call("drawimmediate");
	}
	
	/* use opengl line_smooth command */
	public void linesmooth(int v) {
	    lineSmooth = (v==1) ? true : false;
	}
	
	/* set center position of rainbow */
	public void position(float x, float y, float z) {
		rPosition[0] = x;
		rPosition[1] = y;
		rPosition[2] = z;
	}
	
	/* set main scale of rainbow */
	public void size(float x, float y, float z) {
		rSize[0] = x;
		rSize[1] = y;
		rSize[2] = z;
	}
	
	/* color of the outline, if displayed */
	public void outcolor(float r, float g, float b) {
		outColor[0] = r;
		outColor[1] = g;
		outColor[2] = b;
	}
	
	public void rbase(float x, float y, float z) {
		rBase[0] = x;
		rBase[1] = y;
		rBase[2] = z;
	}
	
	public void rtop(float x, float y, float z) {
		rTop[0] = x;
		rTop[1] = y;
		rTop[2] = z;
	}
	
	/* turn outline on off */
	public void outline(int v) {
		outLine = (v==1) ? true : false;
	}
	
	/* degree of order of bezier curve, the higher the softer */
	public void order(int v) {
		bOrder = (v>0) ? v : 1;
	}
	
	/* resolution of bezier curve */
	public void slices(int v) {
		bSlices = (v>2) ? v : 2;
	}
	
	/* changes the width of the stroke */
	public void scale(float v) {
		scaleRange = (v>0) ? v : 0.01f;
	}
	
	/* toggle, vary stroke width with scaleRange value */
	public void varyscale(int v) {
		varyScale = (v==1) ? true : false;
	}
	
	/* toggle display of bezier points and lines */
	public void showpoints(int v) {
		showPoints = (v==1) ? true : false;
	}
	
	/* toggle wireframe display */
	public void wireframe(int v) {
		wireFrame = (v==1) ? true : false;
	}
	
	/* toggle HSB versus fake rainbow spektrum */
	public void colormode(int v) {
		colorMode = (v==1) ? true : false;
	}
	
	/* return red color component for rainbow band between 0 .. 1 */
	private float rRed(float i) {
		return Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
//		return 1.0f - 1.0f * i;
	}
	
	/* return red color component for rainbow band between 0 .. 1 */
	private float rBlue(float i) {
		return Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
//		return 1.0f * i;
	}
	
	/* return red color component for rainbow band between 0 .. 1 */
	private float rGreen(float i) {
		return Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
//		return 1.0f * i;
	}
	
	private float rainbowSpectrum(float v) {
		if(!colorMode) return v;
		float k = v*2.f;
		float m = (k>1) ? 1 + ((float) Math.sin((k-1.f)*1.570796)) : 1 - (float) Math.cos(k*1.570796);
		return m/2.0f;
	}
	
	/* define number of bands of the rainbow */
	public void bands(int v) {
		rBands = (v>0) ? v : 1;
	}
	
	/* total width of all rainbow bands */
	public void width(float v) {
		rWidth = (v>0.1f) ? v : 0.1f;
	}

	/* set number of bezier points */
	public void pointcount(int v) {
		pointCount = (v>2) ? v : 2;
		if(pointCount>maxPoints) pointCount = maxPoints;
		generateRainbow();
	}
	
	/* toggle morphing */
	public void morph(int v) {
		morphing = (v==1) ? true : false;
	}
	
	/* set morphing speed */
	public void morphspeed(float v) {
		slew = (v>0) ? v : 0.01f;
	}
	
	/* choose design template for rainbow */
	public void mode(int v) {
		rMode = (v>0) ? v : 0;
		generateRainbow();
	}
	
}
