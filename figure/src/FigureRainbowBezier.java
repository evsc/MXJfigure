/*
 * Figure de la terre - mxj objects
 * 
 * make rainbow bezier curves
 * recursive bezier code based on bezier computation based on 
 * http://html5tutorial.com/how-to-draw-n-grade-bezier-curve-with-canvas-api/
 * 
 * animation with addSegment
 * 
 */



import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;


public class FigureRainbowBezier extends MaxObject {

	// render context for opengl objects
	private String context = "foo";
		
	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	private int texture_width = 640;	// TODO: include resize function
	private int texture_height = 480;
		
	// Bezier curve, rainbow backbone
	private int pointCount = 4;
	private int maxPoints = 10;
	private float bp[][];			// bezier points
	private float bc[][];			// bezier colour
	private boolean varyColor = false;
	
	private int bSlices = 20;		// bezier slices
	

	// Bezier display
	private boolean colorMode = true;
	private boolean lineSmooth = true;
	private boolean wireFrame = false;
	private boolean outLine = false;
	private float[] outColor = {0,0,0,1};
	private boolean showPoints = false;
	
	
	// rainbow
	private int rMode = 1;			// rainbow building mode
	private int rBands = 7;
	private float pWidth[];
	private float pScale[];
	private float scaleRange = 1.0f;
	private float widthRange = 0.2f;
	private boolean varyScale = false;
	private boolean varyWidth = false;
	private float[] rPosition = {0,0,0};	// center position of rainbow (arch center on ground)
	private float[] rSize = {1,1,1};		// scale of rainbow
	private float rDir = 1;					// direction. -1 inverts direction of bands

	// position and scale arrays for morphing
	private boolean morphing = false;
	private float slew = 0.1f;
	private float mp[][];				// morph points
	private float mScale[];
	private float mWidth[];
	private float mPosition[] = {0,0,0};
	private float mSize[] = {1,1,1};
	private float mBands = 7;
	private float mScaleRange = 0.1f;
	private float mWidthRange = 0.2f;
	
	// animation
	private boolean autoAnimate = false;
	private float animCounter = 0.f;
	private float animSpeed = 0.01f;
	
	/* instantiating mxj without argument causes a problem */
	public FigureRainbowBezier() {
		// bail method presents instantiation of object and prints error message 
		bail("please provide render context as argument");
	}
	
	
	/* instantiate mxj with render context as argument */
	public FigureRainbowBezier(String c) {
		context = c;			// 	render context
		declareIO(2,1);			/* 	declare number of inlets and outlets, 
 									of DataTypes.ALL */
		
		// assist message for inlets and outlets (for mouse hover)
		setInletAssist(new String[] {"bang to compute and draw", "input settings"});
		setOutletAssist(new String[] {"outputs jit.gl.texture object"});
		
		
		// instantiate jitter objects
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", 1);
		sketch.setAttr("antialias",0);
		sketch.setAttr("glclearcolor", new Atom[]{Atom.newAtom(0.),Atom.newAtom(0.),Atom.newAtom(0.),Atom.newAtom(1.)});
		sketch.setAttr("fsaa", 0);
		sketch.send("automatic", 0);	/* set to not-automatic, to be able to use
										   begin_capture and drawimmediate */

		texture = new JitterObject("jit.gl.texture");
		texture.setAttr("drawto", context);
		texture.setAttr("dim", new Atom[]{Atom.newAtom(texture_width),Atom.newAtom(texture_height)});
		
		// set arrays to max number of pointcount
		int m = maxPoints;
		bp = new float[m][3];
		bc = new float[m][3];
		pScale = new float[m];
		mp = new float[m][3];
		mScale = new float[m];
		pWidth = new float[m];
		mWidth = new float[m];
	}
	
	
	/* defines the points of the rainbow randomly */
	public void generateRainbow() {
		post("generateRainbow() mode:"+rMode);
		
		switch(rMode) {
		case 0:	// random points
			for(int i=0; i<pointCount; i++) {
				bp[i][0] = (float) Math.random()*2-1.f;
				bp[i][1] = (float) Math.random()*2-1.f;
				bp[i][2] = (float) Math.random()*2-1.f;
				pScale[i] = (float) Math.random();
				pWidth[i] = (float) Math.random();
			}
			break;
		case 1: // strict rainbow from 4 points
				pointCount = 4;
				bp[0][0] = bp[1][0] = rPosition[0] - rSize[0];
				bp[2][0] = bp[3][0] = rPosition[0] + rSize[0];
				bp[0][1] = bp[3][1] = rPosition[1];
				bp[1][1] = bp[2][1] = rPosition[1] + rSize[1];
				bp[0][2] = bp[1][2] = bp[2][2] = bp[3][2] = 0;
				for(int i=0; i<pointCount; i++) pScale[i] = 1.0f;
				for(int i=0; i<pointCount; i++) pWidth[i] = 1.0f;
			break;
			
		case 2: // strict rainbow from 3 points
				pointCount = 3;
				bp[0][0] = rPosition[0] - rSize[0];
				bp[2][0] = rPosition[0] + rSize[0];
				bp[1][0] = rPosition[0];
				bp[0][1] = bp[2][1] = rPosition[1];
				bp[1][1] = rPosition[1] + rSize[1];
				bp[0][2] = bp[1][2] = bp[2][2] = 0;
				for(int i=0; i<pointCount; i++) pScale[i] = 1.0f;
				for(int i=0; i<pointCount; i++) pWidth[i] = 1.0f;
			break;
		
		case 3: // advance array of points, one new random point
				for(int i=0; i<pointCount; i++) {
					if(i<pointCount-1) {
						bp[i][0] = bp[i+1][0];
						bp[i][1] = bp[i+1][1];
						bp[i][2] = bp[i+1][2];
					}
					pScale[i] = pScale[i+1];
					pWidth[i] = pWidth[i+1];
				}
				bp[pointCount-1][0] = (float) Math.random()*2-1.f;
				bp[pointCount-1][1] = (float) Math.random()*2-1.f;
				bp[pointCount-1][2] = (float) Math.random()*2-1.f;
				pScale[pointCount-1] = (float) Math.random();
				pWidth[pointCount-1] = (float) Math.random();
			break;
		
			
		}
		
		animCounter = 0.f;
	}
	
	
	/* add new segment to bezier curve, random point */
	public void addSegment() {
		float nx = (float) Math.random()*2-1.f;
		float ny = (float) Math.random()*2-1.f;
		float nz = (float) Math.random()*2-1.f;
		addSegment(nx, ny, nz);
	}
	
	
	/* add new segment to bezier, specific point */
	public void addSegment(float nx, float ny, float nz) {
		
		// need to reduce 
		while(pointCount>maxPoints-3) {
			// strip out first 3 points of array
			for(int i=0; i<pointCount-3; i++) {
				bp[i][0] = bp[i+3][0];
				bp[i][1] = bp[i+3][1];
				bp[i][2] = bp[i+3][2];
				pScale[i] = pScale[i+3];
				pWidth[i] = pWidth[i+3];
			}
			pointCount-=3;
		}
		
		// make sure there are at least 4 points
		if(pointCount<4) {
			pointCount = 4;
			generateRainbow();
		}
		
		// add point 1, compute from second to last point
		int _i = pointCount;
		float[] p_last = { bp[_i-1][0], bp[_i-1][1], bp[_i-1][2]};
		float[] p_last2 = { bp[_i-2][0], bp[_i-2][1], bp[_i-2][2]};
		bp[_i][0] = p_last[0] + (p_last[0]-p_last2[0]);
		bp[_i][1] = p_last[1] + (p_last[1]-p_last2[1]);
		bp[_i][2] = p_last[2] + (p_last[2]-p_last2[2]);
		pScale[_i] = (float) Math.random();
		pWidth[_i] = (float) Math.random();
		
		// add point 2, same as point 1
		_i = pointCount+1;
//		bp[_i][0] = nx + (float) Math.random()-0.5f;
//		bp[_i][1] = ny + (float) Math.random()-0.5f;
//		bp[_i][2] = nz + (float) Math.random()-0.5f;
//		bp[_i][0] = bp[_i-1][0];
//		bp[_i][1] = bp[_i-1][1];
//		bp[_i][2] = bp[_i-1][2];
//		bp[_i][0] = nx + (nx - bp[_i-1][0])/2 + (float) Math.random()*0.1f-0.05f;
//		bp[_i][1] = nx + (ny - bp[_i-1][1])/2 + (float) Math.random()*0.1f-0.05f;
//		bp[_i][2] = nx + (nz - bp[_i-1][2])/2 + (float) Math.random()*0.1f-0.05f;
		float dd = (Math.random() > 0.5f) ? 1.f : -1.f;
		bp[_i][0] = restrict(nx + (bp[_i-1][1] - ny)*dd + (float) Math.random()*0.2f-0.1f);
		bp[_i][1] = restrict(ny -(bp[_i-1][0] - nx)*dd + (float) Math.random()*0.2f-0.1f);
		bp[_i][2] = restrict(nz + (float) Math.random()*0.5f-0.25f);
		pScale[_i] = (float) Math.random();
		pWidth[_i] = (float) Math.random();
		
		// add point 3, specific goal point
		_i = pointCount+2;
		bp[_i][0] = nx;
		bp[_i][1] = ny;
		bp[_i][2] = nz;
		pScale[_i] = (float) Math.random();
		pWidth[_i] = (float) Math.random();
		
		pointCount+=3;

	}
	
	private float restrict(float v) {
		return restrict(v,-1.f,1.f);
	}
	
	private float restrict(float v, float min, float max) {
		if(v<min) return min;
		if(v>max) return max;
		return v;
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
		int _pC = pointCount;
		float _bp[][] = new float[_pC][3];
		float _ps[] = new float[_pC];
		float _pw[] = new float[_pC];
		
		// if morphing, use mx[] as morphed values while px[] represents goal of morphing
		for(int i=0; i<_pC; i++) {
			mp[i][0] = (morphing) ? mp[i][0] += (bp[i][0]-mp[i][0])*slew : bp[i][0];
			mp[i][1] = (morphing) ? mp[i][1] += (bp[i][1]-mp[i][1])*slew : bp[i][1];
			mp[i][2] = (morphing) ? mp[i][2] += (bp[i][2]-mp[i][2])*slew : bp[i][2];
			_bp[i][0] = mp[i][0];
			_bp[i][1] = mp[i][1];
			_bp[i][2] = mp[i][2];
			_ps[i] = (morphing) ? mScale[i] += (pScale[i]-mScale[i])*slew : pScale[i];
			_pw[i] = (morphing) ? mWidth[i] += (pWidth[i]-mWidth[i])*slew : pWidth[i];
		}
		
		// more local variables
		mBands = (morphing) ? mBands += (rBands - mBands)*slew : rBands;
		int _bands = (int) Math.floor(mBands);
		mWidthRange = (morphing) ? mWidthRange += (widthRange - mWidthRange)*slew : widthRange;
		float _sw = mWidthRange;
		mScaleRange = (morphing) ? mScaleRange += (scaleRange - mScaleRange)*slew : scaleRange;
		float _sr = mScaleRange;
		float _add = 1.f / (float) bSlices;			// defines number of slices
		
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
		
		float _d = rDir;	// direction the rainbow expands to
		
		
		sketch.call("reset");
		
		// "line_smooth" 
		if(lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");
		
//		if(wireFrame) sketch.call("glpolygonmode", new Atom[]{Atom.newAtom("front_and_back"),Atom.newAtom("line")});
//		else sketch.call("glpolygonmode", new Atom[]{Atom.newAtom("front_and_back"),Atom.newAtom("fill")});

		

		for(int b=0; b<_bands; b++) {
			
			float index = (_bands==1) ? 0.f : b * 1.0f / (float) (_bands-1); // 0... 1st outer band, 1... inner band
			float cindex = b * 1.0f / (float) _bands;
			
			Atom[] bandColor = rColor(cindex);
		
			// --- - - - - - - - - --- - -start BEZIER - - - -- - - -- - -- - -- -- -- - ---- 
			
	
			sketch.call("glcolor", bandColor);
			sketch.call("gllinewidth", _sr);
			
			int segment = 0;
//			post("segment "+segment+" / pointCount "+_pC+" / pointCount-3 "+(_pC-3));
			while(segment < _pC-3) {
				float _segm[][] = copyArray(_bp,segment,segment+4);
				
				
				if(varyScale) {
//					float _nw = scaleBetween(_ps[segment],_ps[segment+3],t)*_sr;
					float _nw = (0.7f + segment/3.f)*_sr;
					sketch.call("gllinewidth", _nw);	// can only be called outside of glbegin-glend
				}
				sketch.call("glbegin", "line_strip");
				
			
				float st = (segment==0 && autoAnimate) ? animCounter : 0.f;
				float ti = (segment!=0 && segment>=_pC-4 && autoAnimate) ? animCounter : 1.f;
				for(float t=st; t<ti; t+=_add) {
					
	//				if(varyColor) 
					
					float[] _p = P(t, _segm);			// returns array with xyz of main point
					float[] _p2 = P(t+_add, _segm);	// get next point to be able to calculate offset
					float _thew = (varyWidth) ? scaleBetween(_pw[segment],_pw[segment+3],t)*_sw : _sw;
					float[] _off = offsetPoint(_p,_p2,index*_thew*_d);	// calculates offset vector
					
					// combine point with offset vector, depending on which rainbow bow we are drawing right now
					float newx = _x + _p[0]*_sx + _off[0]*_sx;
					float newy = _y + _p[1]*_sy + _off[1]*_sy;
					float newz = _z + _p[2]*_sz + _off[2]*_sz;
					
					sketch.call("glvertex", new Atom[]{Atom.newAtom(newx),Atom.newAtom(newy),Atom.newAtom(newz)});
				}
				
				if(ti>=1.f) {
					// draw last point in curve, with fixed t=1.0 value, to make sure it is precise
					float t=1.0f;
//					if(varyScale) sketch.call("gllinewidth", scaleBetween(_ps[segment],_ps[segment+3],t)*_sr);
					float[] _p = P(t, _segm);	// returns array with xyz of point
					boolean endOfBezier = (segment+4 >= _pC) ? true : false;
					float[] _p2 = endOfBezier ? P(t-_add, _segm) : P(_add, copyArray(_bp,segment+3,segment+7));
					float _thew = (varyWidth) ? scaleBetween(_pw[segment],_pw[segment+3],t)*_sw : _sw;
					float[] _off = offsetPoint(_p,_p2,index*_thew*_d);
					
					// inverse offset on last point of curve, because p2 is the previous instead of the next point in the curve
					float _inv = endOfBezier ? -1.f : 1.f;
					float newx = _x + _p[0]*_sx + _off[0]*_sx*_inv;	
					float newy = _y + _p[1]*_sy + _off[1]*_sy*_inv;
					float newz = _z + _p[2]*_sz + _off[2]*_sz*_inv;
					sketch.call("glvertex", new Atom[]{Atom.newAtom(newx),Atom.newAtom(newy),Atom.newAtom(newz)});
				}
				
				sketch.call("glend");
				segment +=3;
			
			}
			// --- - - - - - - - - --- - - --end BEZIER -- - - -- - -- - -- -- -- - ---- -
		
		}

		// display the control points of the main bezier curve
		if(showPoints) {
			sketch.call("gllinewidth", 1.0f);
			sketch.call("glbegin", "line_strip");
			sketch.call("glcolor", new Atom[]{Atom.newAtom(1.f),Atom.newAtom(1.f),Atom.newAtom(1.f),Atom.newAtom(1.0f)});
			for(int i=0; i<_pC; i++) {
				sketch.call("glvertex", new Atom[]{Atom.newAtom(_x+_bp[i][0]*_sx),Atom.newAtom(_y+_bp[i][1]*_sy),Atom.newAtom(_z+_bp[i][2]*_sz)});
			}
			sketch.call("glend");
			for(int i=0; i<_pC; i++) {
				sketch.call("moveto", new Atom[]{Atom.newAtom(_x+_bp[i][0]*_sx),Atom.newAtom(_y+_bp[i][1]*_sy),Atom.newAtom(_z+_bp[i][2]*_sz)});
				sketch.call("circle", new Atom[]{Atom.newAtom(0.01)});
			}
		}
		// call drawimmediate, to execute drawing of sketch object
		sketch.call("drawimmediate");		
		
		if(autoAnimate) {
			animCounter+=_add;
			if(animCounter>=1.0) {
				animCounter=0.f;
				addSegment();
			}
		}
	}
	
	
	private int tToIndex(float t, int _pC) {
		// t between 0 and 1
		// index = 0 to pointcount - 1
		return (int) (t * (_pC-1));
	}
	
	private float scaleBetween(float v1, float v2, float i) {
//		return (float) Math.random() + 0.1f;
		return v1 + (v2-v1)*i;
	}
	
	private float[][] copyArray(float[][] array0, int from, int to) {
		if(to > array0.length) to = array0.length;
		float[][] r = new float[to-from][3];
		int _i = 0;
		for(int i=from; i<to; i++) {
			r[_i][0] = array0[i][0];
			r[_i][1] = array0[i][1];
			r[_i][2] = array0[i][2];
			_i++;
		}
		return r;
	}
	
	
	/* Bezier : Computes factorial */
	// bezier computation based on http://html5tutorial.com/how-to-draw-n-grade-bezier-curve-with-canvas-api/
	private float fact(float k) {
		if(k==0 || k==1) {
			return 1;
		} else {
			return k* fact(k-1);
		}
	}
	
	/* Bezier : Computes Bernstain
	*@param {Integer} i - the i-th index
	*@param {Integer} n - the total number of points
	*@param {Number} t - the value of parameter t , between 0 and 1
	*/
	private float B(int i, int n, float t) {
		return (float) (fact(n) / (fact(i) * fact(n-i))* Math.pow(t, i) * Math.pow(1-t, n-i));
	}
	
	/* Computes a point's coordinates for a value of t
	*@param {Number} t - a value between o and 1
	*@param {Array} points - an {Array} of [x,y] coordinates. The initial points
	*/
	private float[] P(float t, float[][] points) {
		float[] r = {0,0,0};
		int n = points.length-1;
		for(int i=0; i<=n; i++) {
			r[0] += points[i][0] * B(i,n,t);
			r[1] += points[i][1] * B(i,n,t);
			r[2] += points[i][2] * B(i,n,t);
		}
		return r;
	}

	
	private float[] offsetPoint(float[] p1, float[] p2, float diff) {
		float[] of = {0,0,0};
		float x = p2[0] - p1[0];
		float y = p2[1] - p1[1];
		float z = p2[2] - p1[2];
		float c = (float) Math.sqrt(x*x + y*y + z*z);
		of[0] = -y/c * diff;
		of[1] = x/c * diff;
		of[2] = z/c * diff;
		return of;
	}
	
	
	public void randomizeWidth() {
		for(int i=0; i<pointCount; i++) {
			pScale[i] = (float) Math.random();
			pWidth[i] = (float) Math.random();
		}
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
	
	
	/* turn outline on off */
	public void outline(int v) {
		outLine = (v==1) ? true : false;
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
	
	public void varywidth(int v) {
		varyWidth = (v==1) ? true : false;
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
	
	/* return color of rainbow between 0 ..1 */
	private Atom[] rColor(float i) {
		float r = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
		float g = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
		float b = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
		return new Atom[]{Atom.newAtom(r),Atom.newAtom(g),Atom.newAtom(b)};
	}
	
	/* return red color component for rainbow band between 0 .. 1 */
	private float rRed(float i) {
		return Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
	}
	
	/* return red color component for rainbow band between 0 .. 1 */
	private float rBlue(float i) {
		return Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
	}
	
	/* return red color component for rainbow band between 0 .. 1 */
	private float rGreen(float i) {
		return Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
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
		widthRange = (v>0.1f) ? v : 0.1f;
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
	
	public void animate(int v) {
		autoAnimate = (v==1) ? true : false;
	}
	
		
}
