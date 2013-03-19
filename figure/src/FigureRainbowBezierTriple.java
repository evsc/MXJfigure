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
	private int texture_width = 640;		
	private int texture_height = 240;
		
	// Bezier curve, rainbow backbone
	private int rNum = 1;					// active number of rainbows
	private int rNumMax = 3;				// maximum number of rainbows
	private int pointCount = 4;
	private float bp[][][];					// bezier points
	private float rHeight[];				// height of rainbow arch	
	private int bSlices = 20;				// bezier slices
	

	// Bezier display
	private int colorMode = 1;		// how to calculate rainbow colors

	
	// rainbow
	private int[] rBands = {7,7,7};
	private float[] rWidth = {0.2f,0.2f,0.2f};				// width of all bands together
	private float[] strokeWidth = {0.01f,0.01f,0.01f};		// width of one band
	private float[] rPosition = {0,0,0};	// center position of rainbow (arch center on ground)
	private float[] rSize = {1,1,1};		// scale of rainbow
	private float rDir = 1;					// direction. -1 inverts direction of bands

	// position and scale arrays for morphing
	private int randomMode = 0;				// random building mode
	private boolean morphing = false;
	private float slew = 0.15f;
	private float mp[][][];					// morph points
	private float[] mStrokeWidth = {0.01f,0.01f,0.01f};
	private float mPosition[] = {0,0,0};
	private float mSize[] = {1,1,1};
	private float mWidth[] = { 0.2f,0.2f,0.2f};
	private float[] mBands = {7,7,7};

	
	// display parameters
	private int lineSmooth = 0;
	private int sketchDepthEnable = 1;
	private int sketchAntialias = 0;
	private int sketchFsaa = 0;
	private int sketchBlendEnable = 0;
	private int sketchBlendMode1 = 1;
	private int sketchBlendMode2 = 1;
	
	
	/* instantiating mxj without argument causes a problem */
	public FigureRainbowBezierTriple() {
		// bail method presents instantiation of object and prints error message 
		bail("please provide render context as argument");
	}
	
	
	/* instantiate mxj with render context as argument */
	public FigureRainbowBezierTriple(String c) {
		context = c;			// 	render context
		declareIO(2,2);			/* 	delare number of inlets and outlets, 
 									of DataTypes.ALL */
		
		// assist message for inlets and outlets (for mouse hover)
		setInletAssist(new String[] {"bang to compute and draw", "input settings"});
		setOutletAssist(new String[] {"outputs jit.gl.texture object","connect to thispatcher for gui updating"});
		
		
		// instantiate jitter objects
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", sketchDepthEnable);
		sketch.setAttr("blend_enable", sketchBlendEnable);
		sketch.setAttr("blend_mode", new Atom[] { Atom.newAtom(6), Atom.newAtom(7) });
		sketch.setAttr("antialias", sketchAntialias);
		sketch.setAttr("glclearcolor", new Atom[] { Atom.newAtom(0.), Atom.newAtom(0.),
						Atom.newAtom(1.), Atom.newAtom(1.) });
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
		bp = new float[rNumMax][m][3];
		mp = new float[rNumMax][m][3];
		rHeight = new float[m];
		
		for(int i=0; i<rNumMax; i++) {
			rHeight[i] = 1.f;
			for(int j=0; j<pointCount; j++) {
				for(int k=0; k<3; k++) bp[i][j][k] = 0.f;
			}
		}
		
		generateRainbow();
		outputAllVariables();
	}
	
	
	/* defines the points of the rainbow randomly */
	public void generateRainbow() {
		if(randomMode==1) {
			for(int r=0; r<rNumMax; r++) {
				for(int i=0; i<pointCount; i++) {
					for(int x=0; x<3; x++) {
						bp[r][i][x] = (float) Math.random()*2-1.f;

					}
				}
			}
		} else {
			pA(-2.f,0,-3,-2.1f,0,-3,-1.9f,0.f,-3);
			pB(0.3f,-0.5f,0.4f,-0.2f,-0.3f,1.1f,0.5f,-0.2f,0.7f);
		}
		outputBezierVariables();
	}
	
	private void setBezierPoints() {
		for(int r=0; r<rNumMax; r++) {
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
		if(debug) post("end_capture");
		texture.call("end_capture");			// end capturing
		if(debug) post("draw");
		texture.call("draw");					// to output texture? 
		if(debug) post("jit_gl_texture");
		outlet(0,"jit_gl_texture",texture.getAttr("name"));		// output texture
	}

	
	/* draw rainbow to jitter sketch object */
	public void draw() {


		// create local variables, to avoid conflict when life-updating variables while rendering
		float _bp[][][] = new float[rNumMax][pointCount][3];
		
		sketch.call("reset");	// start drawing by resetting sketch object
		
		if (lineSmooth>0) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");

		for(int r=0; r<rNum; r++) {

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
			mBands[r] = (morphing) ? mBands[r] += (rBands[r] - mBands[r])*slew : rBands[r];
			int _bands = (int) Math.floor(mBands[r]);
			mWidth[r] = (morphing) ? mWidth[r] += (rWidth[r] - mWidth[r])*slew : rWidth[r];
			float _w = mWidth[r];
			mStrokeWidth[r] = (morphing) ? mStrokeWidth[r] += (strokeWidth[r] - mStrokeWidth[r])*slew : strokeWidth[r];
			float _scw = mStrokeWidth[r];
			if(_scw > _w/ ((_bands-1)*2.f)) _scw = _w/ ((_bands-1)*2.f);
			float _add = 1.f / (float) (bSlices);			// defines number of slices

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


			/* calculate all bezier points and their offset vectors beforehand
			 * as they are used multiple times in the following loops */
			float[][] bezierPoint = new float[bSlices+1][3];
			float[][] offsetVector = new float[bSlices+1][3];
			// calculate all points along the bezier curve
			for(int _t=0; _t<=bSlices; _t++) {
				float t = _t*_add;
				bezierPoint[_t] = P(t, _bp[r]);	
			}
			// calculate vector perpendicular to the line connecting the previous and next bezier point
			for(int _t=0; _t<=bSlices; _t++) {
				if(_t==0) offsetVector[_t] = perpendicularVector(bezierPoint[_t],bezierPoint[_t+1]);
				else if(_t==bSlices) offsetVector[_t] = perpendicularVector(bezierPoint[_t-1],bezierPoint[_t]);
				else offsetVector[_t] = perpendicularVector(bezierPoint[_t-1],bezierPoint[_t+1]);
			}

			// calculate offsetPoints of bands before drawing.
			for(int b=0; b<_bands; b++) {

				float index = b * 1.0f / (float) (_bands-1); // 0... 1st outer band, 1... inner band
				float cindex =  b * 1.0f / (float) _bands;

				
				
				float[][][] bezierOffset = new float[bSlices+1][2][3];	// hold offset points on both sides of curve
				for(int _t=0; _t<=bSlices; _t++) {
					float[] _offOutside = scaleVector(offsetVector[_t], index*_w*_d + _scw);
					float[] _offInside = scaleVector(offsetVector[_t], index*_w*_d - _scw);
					
					bezierOffset[_t][0][0] = _x + bezierPoint[_t][0]*_sx + _offOutside[0]*_sx;
					bezierOffset[_t][1][0] = _x + bezierPoint[_t][0]*_sx + _offInside[0]*_sx;
					bezierOffset[_t][0][1] = _y + bezierPoint[_t][1]*_sy + _offOutside[1]*_sy;
					bezierOffset[_t][1][1] = _y + bezierPoint[_t][1]*_sy + _offInside[1]*_sy;
					bezierOffset[_t][0][2] = _z + bezierPoint[_t][2]*_sz + _offOutside[2]*_sz;
					bezierOffset[_t][1][2] = _z + bezierPoint[_t][2]*_sz + _offInside[2]*_sz;
				}
				
				

				// --- - - - - - - - - --- - -start BEZIER - - - -- - - -- - -- - -- -- -- - ---- 

				Atom[] bandColor = rColor(cindex);
				sketch.call("glcolor", bandColor);

				sketch.call("glbegin", "tri_strip");

				for(int i=0; i<=bSlices; i++) {
					sketch.call("glvertex", new Atom[]{	Atom.newAtom(bezierOffset[i][0][0]),
							Atom.newAtom(bezierOffset[i][0][1]),
							Atom.newAtom(bezierOffset[i][0][2])});
					sketch.call("glvertex", new Atom[]{	Atom.newAtom(bezierOffset[i][1][0]),
							Atom.newAtom(bezierOffset[i][1][1]),
							Atom.newAtom(bezierOffset[i][1][2])});
				}
				
				sketch.call("glend");

				// --- - - - - - - - - --- - - --end BEZIER -- - - -- - -- - -- -- -- - ---- -

			}

		}
		
		// call drawimmediate, to execute drawing of sketch object
//		try {
			if(debug) post("drawimmediate");
			sketch.call("drawimmediate");	
//		} catch(Exception e) {
//			if(debug) post("drawimmediate error: "+e);
//		}	

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

	
	private float[] perpendicularVector(float[] point1, float[] point2) {
		
		float[] ofv = {0,0,0};
		
		float x = point2[0] - point1[0];
		float y = point2[1] - point1[1];
		float z = point2[2] - point1[2];
		
		float c = (float) Math.sqrt(x*x + y*y + z*z);
		ofv[0] = -y/c;
		ofv[1] = x/c;
		ofv[2] = z/c;
		
		return ofv;
	}
	
	
	private float[] scaleVector(float[] v, float s) {
		float[] newVector = { v[0]*s, v[1]*s, v[2]*s };
		return newVector;
	}
	
	
	/* return color of rainbow between 0 ..1 */
	private Atom[] rColor(float i) {
		float r = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
		float g = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
		float b = Color.getHSBColor(rainbowSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
		return new Atom[]{Atom.newAtom(r),Atom.newAtom(g),Atom.newAtom(b)};
	}
	
	private float rainbowSpectrum(float v) {
		if(colorMode==0) return v;
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
	
	public void resize(int x, int y) {
		texture_width = x;
		texture_height = y;
		texture.setAttr("dim",new Atom[] { Atom.newAtom(texture_width), Atom.newAtom(texture_height) });
		if(debug) post("resized to "+x + " "+y);
	}
	
	public void depthEnable(int v) {
		sketchDepthEnable = (v == 1) ? 1 : 0;
		sketch.setAttr("depth_enable", sketchDepthEnable);
	}

	public void blendEnable(int v) {
		sketchBlendEnable = (v == 1) ? 1 : 0;
		sketch.setAttr("blend_enable", sketchBlendEnable);
	}

	public void blendMode(int b1, int b2) {
		sketchBlendMode1 = b1;
		sketchBlendMode2 = b2;
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
		lineSmooth = (v == 1) ? 1 : 0;
	}
	
	
	/* === === === === === === === rainbow === === === === === === === === */
	
	
	public void num(int v) {
		rNum = (v<rNumMax) ? (v>0) ? v : 1 : rNumMax;
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
	
	/* set the start points of the rainbow arches */
	public void pA(float x, float y, float z) {
		pA(x,y,z,x,y,z,x,y,z);
	}
	
	/* set the start points of the rainbow arches */
	public void pA(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
		bp[0][0][0] = x1;
		bp[0][0][1] = y1;
		bp[0][0][2] = z1;
		bp[1][0][0] = x2;
		bp[1][0][1] = y2;
		bp[1][0][2] = z2;
		bp[2][0][0] = x3;
		bp[2][0][1] = y3;
		bp[2][0][2] = z3;
		setBezierPoints();
	}

	/* set the end points of the rainbow arches */
	public void pB(float x, float y, float z) {
		pB(x,y,z,x,y,z,x,y,z);
	}
	
	/* set the end points of the rainbow arches */
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
		strokewidth(v,v,v);
	}
	
	public void strokewidth(float v1, float v2, float v3) {
		strokeWidth[0] = (v1>0) ? v1 : 0.001f;
		strokeWidth[1] = (v2>0) ? v2 : 0.001f;
		strokeWidth[2] = (v3>0) ? v3 : 0.001f;
	}
	
	
	/* toggle HSB versus fake rainbow spektrum */
	public void colormode(int v) {
		colorMode = (v==1) ? 1 : 0;
	}

	
	/* define number of bands of the rainbow */
	public void bands(int v) {
		bands(v,v,v);
	}
	
	public void bands(int v1, int v2, int v3) {
		rBands[0] = (v1>0) ? v1 : 1;
		rBands[1] = (v2>0) ? v2 : 1;
		rBands[2] = (v3>0) ? v3 : 1;
	}
	
	/* total width of all rainbow bands */
	public void width(float v) {
		width(v,v,v);
	}
	
	public void width(float v1, float v2, float v3) {
		rWidth[0] = v1;
		rWidth[1] = v2;
		rWidth[2] = v3;
	}
	
	/* set the height of the rainbow arches */
	public void height(float v) {
		height(v,v,v);
	}
	
	/* set the height of the rainbow arches */
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
		randomMode = (v==1) ? 1 : 0;
		generateRainbow();
	}
	
	
	
	/* === === === === === === === GUI in max patch === === === === === === === === */
	
	public void outputVariables() {
		outputAllVariables();
		outputBezierVariables();
	}
	
	private void outputAllVariables() {

		// display parameters
		outlet(1,"script send gui_depthenable set "+sketchDepthEnable);
		outlet(1,"script send gui_antialias set "+sketchAntialias);
		outlet(1,"script send gui_fsaa set "+sketchFsaa);
		outlet(1,"script send gui_linesmooth set "+lineSmooth);
		outlet(1,"script send gui_blendenable set "+sketchBlendEnable);
		outlet(1,"script send gui_blendmode1 set "+sketchBlendMode1);
		outlet(1,"script send gui_blendmode2 set "+sketchBlendMode2);
		
		// rainbow
		outlet(1,"script send gui_colormode set "+colorMode);
		outlet(1,"script send gui_sizex set "+rSize[0]);
		outlet(1,"script send gui_sizey set "+rSize[1]);
		outlet(1,"script send gui_sizez set "+rSize[2]);
		outlet(1,"script send gui_positionx set "+rPosition[0]);
		outlet(1,"script send gui_positiony set "+rPosition[1]);
		outlet(1,"script send gui_positionz set "+rPosition[2]);
		
		
		outlet(1,"script send gui_slices set "+bSlices);
		
		outlet(1,"script send gui_bands0 set "+rBands[0]);
		outlet(1,"script send gui_bands1 set "+rBands[0]);
		outlet(1,"script send gui_bands2 set "+rBands[1]);
		outlet(1,"script send gui_bands3 set "+rBands[2]);
		outlet(1,"script send gui_bands set bands "+rBands[0]+" "+rBands[1]+" "+rBands[2]);
		
		outlet(1,"script send gui_width0 set "+rWidth[0]);
		outlet(1,"script send gui_width1 set "+rWidth[0]);
		outlet(1,"script send gui_width2 set "+rWidth[1]);
		outlet(1,"script send gui_width3 set "+rWidth[2]);
		outlet(1,"script send gui_width set width "+rWidth[0]+" "+rWidth[1]+" "+rWidth[2]);
		
		outlet(1,"script send gui_strokewidth0 set "+strokeWidth[0]);
		outlet(1,"script send gui_strokewidth1 set "+strokeWidth[0]);
		outlet(1,"script send gui_strokewidth2 set "+strokeWidth[1]);
		outlet(1,"script send gui_strokewidth3 set "+strokeWidth[2]);
		outlet(1,"script send gui_strokewidth set strokewidth "+strokeWidth[0]+" "+strokeWidth[1]+" "+strokeWidth[2]);
		
		outlet(1,"script send gui_height0 set "+rHeight[0]);
		outlet(1,"script send gui_height1 set "+rHeight[0]);
		outlet(1,"script send gui_height2 set "+rHeight[1]);
		outlet(1,"script send gui_height3 set "+rHeight[2]);
		outlet(1,"script send gui_height set height "+rHeight[0]+" "+rHeight[1]+" "+rHeight[2]);
		
		outlet(1,"script send gui_morphspeed set "+slew);
		outlet(1,"script send gui_morph set "+ (morphing ? 1 : 0));
		outlet(1,"script send gui_mode set "+randomMode);
		
	}
	
	private void outputBezierVariables() {
		outlet(1,"script send gui_pa1x set "+bp[0][0][0]);
		outlet(1,"script send gui_pa1y set "+bp[0][0][1]);
		outlet(1,"script send gui_pa1z set "+bp[0][0][2]);
		outlet(1,"script send gui_pa2x set "+bp[1][0][0]);
		outlet(1,"script send gui_pa2y set "+bp[1][0][1]);
		outlet(1,"script send gui_pa2z set "+bp[1][0][2]);
		outlet(1,"script send gui_pa3x set "+bp[2][0][0]);
		outlet(1,"script send gui_pa3y set "+bp[2][0][1]);
		outlet(1,"script send gui_pa3z set "+bp[2][0][2]);
		
		outlet(1,"script send gui_pa set pA "
				+bp[0][0][0]+" "+bp[0][0][1]+" "+bp[0][0][2]+" "
				+bp[1][0][0]+" "+bp[1][0][1]+" "+bp[1][0][2]+" "
				+bp[2][0][0]+" "+bp[2][0][1]+" "+bp[2][0][2]+" ");
		
		
		outlet(1,"script send gui_pb1x set "+bp[0][pointCount - 1][0]);
		outlet(1,"script send gui_pb1y set "+bp[0][pointCount - 1][1]);
		outlet(1,"script send gui_pb1z set "+bp[0][pointCount - 1][2]);
		outlet(1,"script send gui_pb2x set "+bp[1][pointCount - 1][0]);
		outlet(1,"script send gui_pb2y set "+bp[1][pointCount - 1][1]);
		outlet(1,"script send gui_pb2z set "+bp[1][pointCount - 1][2]);
		outlet(1,"script send gui_pb3x set "+bp[2][pointCount - 1][0]);
		outlet(1,"script send gui_pb3y set "+bp[2][pointCount - 1][1]);
		outlet(1,"script send gui_pb3z set "+bp[2][pointCount - 1][2]);
		
		outlet(1,"script send gui_pb set pB "
				+bp[0][pointCount - 1][0]+" "+bp[0][pointCount - 1][1]+" "+bp[0][pointCount - 1][2]+" "
				+bp[1][pointCount - 1][0]+" "+bp[1][pointCount - 1][1]+" "+bp[1][pointCount - 1][2]+" "
				+bp[2][pointCount - 1][0]+" "+bp[2][pointCount - 1][1]+" "+bp[2][pointCount - 1][2]+" ");
	}
	
		
}
