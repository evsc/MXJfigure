/*
 * Figure de la terre - mxj objects
 * 
 * make rainbow bezier curves
 * recursive bezier code based on bezier computation based on 
 * http://html5tutorial.com/how-to-draw-n-grade-bezier-curve-with-canvas-api/
 * 
 * triple rainbow, no animation
 * 
 */



import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;


public class FigureRainbowBezierTriple extends MaxObject {

	// render context for opengl objects
	private String context = "foo";
	private boolean debug = false;
		
	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	private int texture_width = 640;	// TODO: include resize function
	private int texture_height = 480;
		
	// Bezier curve, rainbow backbone
	private int rNo = 3;			// number of rainbows
	private int pointCount = 4;
	private float bp[][][];			// bezier points
	private float rHeight[];		// height of rainbow arch
	private float pScale[];
	private boolean varyScale = false;
	private float scaleRange = 0.1f;
	private int bSlices = 20;		// bezier slices
	

	// Bezier display
	private boolean colorMode = true;		// how to calculate rainbow colors

	
	// rainbow
	private int rMode = 1;					// rainbow building mode
	private int rBands = 7;
	private float rWidth = 0.2f;
	private float[] rPosition = {0,0,0};	// center position of rainbow (arch center on ground)
	private float[] rSize = {1,1,1};		// scale of rainbow
	private float rDir = 1;					// direction. -1 inverts direction of bands

	// position and scale arrays for morphing
	private boolean morphing = false;
	private float slew = 0.1f;
	private float mp[][][];				// morph points
	private float mScale[];
	private float mPosition[] = {0,0,0};
	private float mSize[] = {1,1,1};
	private float mWidth = 0.2f;
	private float mBands = 7;
	private float mScaleRange = 0.1f;
	
	
	// display parameters
	private boolean lineSmooth = false;
	private int sketchDepthEnable = 0;
	private int sketchAntialias = 0;
	private int sketchFsaa = 0;
	private int sketchBlendEnable = 0;
	
	
	
	/* instantiating mxj without argument causes a problem */
	public FigureRainbowBezierTriple() {
		// bail method presents instantiation of object and prints error message 
		bail("please provide render context as argument");
	}
	
	
	/* instantiate mxj with render context as argument */
	public FigureRainbowBezierTriple(String c) {
		context = c;			// 	render context
		declareIO(2,1);			/* 	delare number of inlets and outlets, 
 									of DataTypes.ALL */
		
		// assist message for inlets and outlets (for mouse hover)
		setInletAssist(new String[] {"bang to compute and draw", "input settings"});
		setOutletAssist(new String[] {"outputs jit.gl.texture object"});
		
		
		// instantiate jitter objects
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", sketchDepthEnable);
		sketch.setAttr("blend_enable", sketchBlendEnable);
		sketch.setAttr("blend_mode", new Atom[] { Atom.newAtom(6), Atom.newAtom(7) });
		sketch.setAttr("antialias", sketchAntialias);
		sketch.setAttr("glclearcolor", new Atom[] { Atom.newAtom(0.), Atom.newAtom(0.),
						Atom.newAtom(0.), Atom.newAtom(1.) });
		sketch.setAttr("fsaa", sketchFsaa);
		sketch.send("automatic", 0); /*
									 * set to not-automatic, to be able to use
									 * begin_capture and drawimmediate for
									 * capturing jit.gl.sketch as texture
									 */

		texture = new JitterObject("jit.gl.texture");
		texture.setAttr("drawto", context);
		texture.setAttr("dim", new Atom[]{Atom.newAtom(texture_width),Atom.newAtom(texture_height)});
		
		// set arrays to max number of pointcount
		int m = pointCount;
		bp = new float[rNo][m][3];
		pScale = new float[m];
		mp = new float[rNo][m][3];
		mScale = new float[m];
		rHeight = new float[m];
		
		for(int i=0; i<rNo; i++) {
			rHeight[i] = 1.f;
			pScale[i] = 1.f;
			for(int j=0; j<pointCount; j++) {
				for(int k=0; k<3; k++) bp[i][j][k] = 0.f;
			}
		}
		
		pA(0,0,-3);
		pB(-1,0,1,0,0,1,1,0,1);
	}
	
	
	/* defines the points of the rainbow randomly */
	public void generateRainbow() {
		post("generateRainbow() mode:"+rMode);
		
		switch(rMode) {
		case 0:	// random points
			for(int r=0; r<rNo; r++) {
				for(int i=0; i<pointCount; i++) {
					for(int x=0; x<3; x++) {
						bp[r][i][x] = (float) Math.random()*2-1.f;

					}
					pScale[i] = (float) Math.random();
				}
			}
			break;
		case 1: // strict rainbow from 4 points
			for(int r=0; r<rNo; r++) {
				bp[r][0][0] = bp[r][1][0] = -rSize[0];
				bp[r][2][0] = bp[r][3][0] = rSize[0];
				bp[r][0][1] = bp[r][3][1] = 0;
				bp[r][1][1] = bp[r][2][1] = rSize[1];
				bp[r][0][2] = bp[r][1][2] = bp[r][2][2] = bp[r][3][2] = 0;
			}
			for(int i=0; i<pointCount; i++) pScale[i] = 1.0f;
			break;
		}
	}
	
	private void setBezierPoints() {
		for(int r=0; r<rNo; r++) {
			bp[r][1][0] = bp[r][0][0];
			bp[r][1][1] = bp[r][0][1] + rHeight[r];
			bp[r][1][2] = bp[r][0][2];
			
			bp[r][2][0] = bp[r][3][0];
			bp[r][2][1] = bp[r][3][1] + rHeight[r];
			bp[r][2][2] = bp[r][3][2];
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
		float _bp[][][] = new float[rNo][pointCount][3];
		
		sketch.call("reset");	// start drawing by resetting sketch object
		
		if (lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");

		for(int r=0; r<rNo; r++) {

			// if morphing, use mx[] as morphed values while px[] represents goal of morphing
			for(int i=0; i<pointCount; i++) {
				mp[r][i][0] = (morphing) ? mp[r][i][0] += (bp[r][i][0]-mp[r][i][0])*slew : bp[r][i][0];
				mp[r][i][1] = (morphing) ? mp[r][i][1] += (bp[r][i][1]-mp[r][i][1])*slew : bp[r][i][1];
				mp[r][i][2] = (morphing) ? mp[r][i][2] += (bp[r][i][2]-mp[r][i][2])*slew : bp[r][i][2];
				_bp[r][i][0] = mp[r][i][0];
				_bp[r][i][1] = mp[r][i][1];
				_bp[r][i][2] = mp[r][i][2];
			}

			// more local variables
			mBands = (morphing) ? mBands += (rBands - mBands)*slew : rBands;
			int _bands = (int) Math.floor(mBands);
			mWidth = (morphing) ? mWidth += (rWidth - mWidth)*slew : rWidth;
			float _w = mWidth;
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

			//		if(wireFrame) sketch.call("glpolygonmode", new Atom[]{Atom.newAtom("front_and_back"),Atom.newAtom("line")});
			//		else sketch.call("glpolygonmode", new Atom[]{Atom.newAtom("front_and_back"),Atom.newAtom("fill")});



			for(int b=0; b<_bands; b++) {

				float index = b * 1.0f / (float) (_bands-1); // 0... 1st outer band, 1... inner band
				float cindex =  b * 1.0f / (float) _bands;

				Atom[] bandColor = rColor(cindex);

				// --- - - - - - - - - --- - -start BEZIER - - - -- - - -- - -- - -- -- -- - ---- 


				sketch.call("glcolor", bandColor);
				sketch.call("gllinewidth", pScale[0]*_sr);

				sketch.call("glbegin", "line_strip");

				for(float t=0.0f; t<1.0f; t+=_add) {

					float[] _p = P(t, _bp[r]);			// returns array with xyz of main point
					float[] _p2 = P(t+_add, _bp[r]);	// get next point to be able to calculate offset
					float[] _off = offsetPoint(_p,_p2,index*_w*_d);	// calculates offset vector

					// combine point with offset vector, depending on which rainbow bow we are drawing right now
					float newx = _x + _p[0]*_sx + _off[0]*_sx;
					float newy = _y + _p[1]*_sy + _off[1]*_sy;
					float newz = _z + _p[2]*_sz + _off[2]*_sz;

					sketch.call("glvertex", new Atom[]{Atom.newAtom(newx),Atom.newAtom(newy),Atom.newAtom(newz)});
				}

				// draw last point in curve, with fixed t=1.0 value, to make sure it is precise
				float t=1.0f;

				float[] _p = P(t, _bp[r]);	// returns array with xyz of point
				float[] _p2 = P(t-_add, _bp[r]);
				float[] _off = offsetPoint(_p,_p2,index*_w*_d);
				// inverse offset on last point of curve, because p2 is the previous instead of the next point in the curve
				float newx = _x + _p[0]*_sx - _off[0]*_sx;	
				float newy = _y + _p[1]*_sy - _off[1]*_sy;
				float newz = _z + _p[2]*_sz - _off[2]*_sz;
				sketch.call("glvertex", new Atom[]{Atom.newAtom(newx),Atom.newAtom(newy),Atom.newAtom(newz)});

				sketch.call("glend");

				// --- - - - - - - - - --- - - --end BEZIER -- - - -- - -- - -- -- -- - ---- -

			}

		}
		
		// call drawimmediate, to execute drawing of sketch object
		sketch.call("drawimmediate");		

	}
	
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === == helpful little functions == === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	private int tToIndex(float t) {
		// t between 0 and 1
		// index = 0 to pointcount - 1
		return (int) (t * (pointCount-1));
	}
	
	
	
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === ===  bezier functions  === === === === === === === 
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
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
	
	
	/* return color of rainbow between 0 ..1 */
	private Atom[] rColor(float i) {
		float r = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
		float g = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
		float b = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
		return new Atom[]{Atom.newAtom(r),Atom.newAtom(g),Atom.newAtom(b)};
	}
	
	private float rainbowSpectrum(float v) {
		if(!colorMode) return v;
		float k = v*2.f;
		float m = (k>1) ? 1 + ((float) Math.sin((k-1.f)*1.570796)) : 1 - (float) Math.cos(k*1.570796);
		return m/2.0f;
	}
	
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === === === set parameters === === === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	/* === === === === === === === jit.gl.sketch === === === === === === === === */
	
	public void depthEnable(int v) {
		sketchDepthEnable = (v == 1) ? 1 : 0;
		sketch.setAttr("depth_enable", sketchDepthEnable);
	}

	public void blendEnable(int v) {
		sketchBlendEnable = (v == 1) ? 1 : 0;
		sketch.setAttr("blend_enable", sketchBlendEnable);
	}

	public void blendMode(int b1, int b2) {
		sketch.setAttr("blend_mode",
				new Atom[] { Atom.newAtom(b1), Atom.newAtom(b2) });
	}

	public void antialias(int v) {
		sketchAntialias = (v == 1) ? 1 : 0;
		sketch.setAttr("antialias", sketchAntialias);
	}

	public void clearColor(float r, float g, float b, float a) {
		sketch.setAttr("glclearcolor",
				new Atom[] { Atom.newAtom(r), Atom.newAtom(g), Atom.newAtom(b),
						Atom.newAtom(a) });
	}

	public void fsaa(int v) {
		sketchFsaa = (v == 1) ? 1 : 0;
		sketch.setAttr("fsaa", sketchFsaa);
	}

	public void linesmooth(int v) {
		lineSmooth = (v == 1) ? true : false;
	}
	
	
	/* === === === === === === === rainbow === === === === === === === === */
	
	
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
	
	
	public void pA(float x, float y, float z) {
		pA(x,y,z,x,y,z,x,y,z);
	}
	
	public void pA(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
		bp[0][0][0] = x1;
		bp[0][0][1] = y1;
		bp[0][0][2] = z1;
		bp[1][0][0] = x2;
		bp[1][0][1] = y2;
		bp[1][0][2] = z2;
		bp[2][0][0] = x3;
		bp[2][0][1] = y2;
		bp[2][0][2] = z3;
		setBezierPoints();
	}

	public void pB(float x, float y, float z) {
		pB(x,y,z,x,y,z,x,y,z);
	}
	
	public void pB(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
		bp[0][pointCount - 1][0] = x1;
		bp[0][pointCount - 1][1] = y1;
		bp[0][pointCount - 1][2] = z1;
		bp[1][pointCount - 1][0] = x2;
		bp[1][pointCount - 1][1] = y2;
		bp[1][pointCount - 1][2] = z2;
		bp[2][pointCount - 1][0] = x3;
		bp[2][pointCount - 1][1] = y3;
		bp[2][pointCount - 1][2] = z3;
		setBezierPoints();
	}
	

	/* resolution of bezier curve */
	public void slices(int v) {
		bSlices = (v>2) ? v : 2;
	}
	
	/* changes the width of the stroke */
	public void strokewidth(float v) {
		scaleRange = (v>0) ? v : 0.01f;
	}
	
	/* toggle, vary stroke width with scaleRange value */
	public void varyscale(int v) {
		varyScale = (v==1) ? true : false;
	}
	
	/* toggle HSB versus fake rainbow spektrum */
	public void colormode(int v) {
		colorMode = (v==1) ? true : false;
	}

	
	/* define number of bands of the rainbow */
	public void bands(int v) {
		rBands = (v>0) ? v : 1;
	}
	
	/* total width of all rainbow bands */
	public void width(float v) {
		rWidth = (v>0.1f) ? v : 0.1f;
	}
	
	public void height(float v) {
		height(v,v,v);
	}
	
	public void height(float v1, float v2, float v3) {
		rHeight[0] = v1;
		rHeight[1] = v2;
		rHeight[2] = v3;
		setBezierPoints();
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
