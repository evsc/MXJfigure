/*
 * Figure de la terre - mxj objects
 * 
 * paisley wallpaper
 * simple grid, bezier element
 * 
 */

import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;



public class FigurePaisley extends MaxObject {
	
	// render context for jitter opengl objects
	private String context = "foo";
	private boolean debug = false;

	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	private int texture_width = 640; 		
	private int texture_height = 240;
	
	
	
	// grid
	private float[] pGrid = { 0.3f, 0.3f };
	private float[] pPos = {0,0,0};	// center position
	private float[] pSize = {2,2,2};		// size of the grid?? 
	private float mirrorDistance = 0.5f;
	
	// 
	private int pNumMax = 24;	// maximum number of paisley patterns
	
	
	// pattern parameter
	private float[][][] masterPattern;		// the master pattern template, 
											// every displayed pattern will be a slight variation of this one
	
	private int eNum = 4;					// number of elements that constitute one paisley pattern
	private int eNumMax = 5;
	private int[] eAnchor;					/* special anchoring of element 
	 										*  0 = no effect
	 										*  1 = keep first point fixed
	 										*  2 = first point is relative to startpoint of previous element
	 										*  3 = first point is relative to endpoint of previous element
	 										*/
											
	private int bNum = 6;					// bezier curves made from maximum n points
	private int[] ePoints;
	private int bSlices = 20;				// bezier slices
	
	private float[][][][] patternNoise;		// noise deviation from the master pattern
	private float bNoiseWeight = 0.05f;		// weight of bezier noise deviation
	
	private float[] masterWidth;			// the curve width of the master pattern
	private float[][] widthNoise;			// noise deviation from master pattern width
	private float wNoiseWeight = 0.1f;		// weight of width noise deviation
	
	
	// morphing
	private boolean morphing = false;		//
	private boolean animation = false;		// random noise changes every frame
	private float noiseFactor = 1.f;		// how likely a noisevalue will be changed
	private float returnFactor = 10.f;		// how likely the noisevalue will be set to 0
	private float slew = 0.1f; 				// morphing speed
	private float[] mPos = { 0, 0, 0 };
	private float[] mSize = { 2, 2, 2 };
	private float[] mGrid = { 0.4f, 0.4f };
	private float mMirrorDistance = 0.1f;
	private float sizeMult = 0.1f;			// scale down size, as 100% would cover whole screen
	private float[][][][] mPatternNoise;	// morphing the noise deviation
	private float[][] mWidthNoise;			// morphing noise width deviation
	
	
	// display parameters
	private boolean lineSmooth = false;
	private int sketchDepthEnable = 0;
	private int sketchAntialias = 0;
	private int sketchFsaa = 0;
	private int sketchBlendEnable = 0;
	private int sketchBlendMode1 = 1;
	private int sketchBlendMode2 = 1;
	
	
	
	/*
	 * === === === === === === === main functions === === === === === === ===
	 */

	
	/* instantiating mxj without argument */
	public FigurePaisley() {
		bail("please provide render context as argument");
	}


	
	/* instantiate mxj with render context as argument */
	public FigurePaisley(String rc) {

		context = rc; 			// render context
		declareIO(2, 2); 		// declare 2 inlets, 1 outlet of DataTypes.ALL

		// assist message for inlets and outlets (mouse hover)
		setInletAssist(new String[] { "bang to compute and draw", "input settings" });
		setOutletAssist(new String[] { "outputs jit.gl.texture object","connect to thispatcher for gui updating" });

		// instantiate Jitter sketch object
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", sketchDepthEnable);
		sketch.setAttr("blend_enable", sketchBlendEnable);
		sketch.setAttr("blend_mode", new Atom[] { Atom.newAtom(6), Atom.newAtom(7) });
		sketch.setAttr("antialias", sketchAntialias);
		sketch.setAttr("glclearcolor",
				new Atom[] { Atom.newAtom(0.), Atom.newAtom(1.), Atom.newAtom(0.), Atom.newAtom(1.) });
		sketch.setAttr("fsaa", sketchFsaa);
		sketch.send("automatic", 0); /*
									 * set to not-automatic, to be able to use
									 * begin_capture and drawimmediate for
									 * capturing jit.gl.sketch as texture
									 */

		// instantiate Jitter texture object
		texture = new JitterObject("jit.gl.texture");
		texture.setAttr("drawto", context);
		texture.setAttr("dim",new Atom[] { Atom.newAtom(texture_width), Atom.newAtom(texture_height) });
		
		
		masterPattern = new float[eNumMax][bNum][3];	// even though we only use xy, keep xyz for consistency
		patternNoise = new float[pNumMax][eNumMax][bNum][3];
		mPatternNoise = new float[pNumMax][eNumMax][bNum][3];
		masterWidth = new float[eNumMax];
		widthNoise = new float[pNumMax][eNumMax];
		mWidthNoise = new float[pNumMax][eNumMax];
		ePoints = new int[eNumMax];
		eAnchor = new int[eNumMax];
		
		for(int e=0; e<eNumMax; e++ ) {
			masterWidth[e] = 0.05f;
			ePoints[e] = bNum;
			eAnchor[e] = 0;
			for(int b=0; b<bNum; b++) {
			}
		}
		for(int p=0; p<pNumMax; p++) {
			for(int e=0; e<eNumMax; e++ ) {
				for(int b=0; b<bNum; b++) {
					mPatternNoise[p][e][b] = new float[] { 0.f, 0.f, 0.f };
				}
				mWidthNoise[p][e] = 0.f;
			}
		}
		
		createMasterPaisley();
		randomize(0.1f,0.f);
	}



	/*
	 * draws and captures the aurora drawing to sketch object, and outputs
	 * sketch as jitter texture object
	 */
	public void bang() {
//		if(debug) post("begin_capture");
		texture.call("begin_capture"); // begin capturing	
//		if(debug) post("draw");
		draw(); // draw aurora to sketch
//		if(debug) post("end_capture");
		texture.call("end_capture"); // end capturing
//		if(debug) post("draw");
		texture.call("draw"); // to output texture?
		outlet(0, "jit_gl_texture", texture.getAttr("name")); // output texture
	}
	
	
	/* draw aurora to jitter sketch object */
	public void draw() {
		
		if(animation) randomize(noiseFactor/100.f, returnFactor/100.f);

//		if(debug) post("reset");
		sketch.call("reset");
		
//		if(debug) post("line_smooth");
		if (lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");
		
		// create local variables, to avoid conflict when life-updating
		// variables while rendering
		
		float _size[] = new float[2];		// scale from origin
		float _pos[] = new float[2];		// origin position
		float _grid[] = new float[2];
		
		for (int i = 0; i < 2; i++) {		// only need x and y values
			mSize[i] = (morphing) ? mSize[i] += (pSize[i] - mSize[i]) * slew : pSize[i];
			_size[i] = mSize[i] * sizeMult;
			mPos[i] = (morphing) ? mPos[i] += (pPos[i] - mPos[i]) * slew : pPos[i];
			_pos[i] = mPos[i];
			mGrid[i] = (morphing) ? mGrid[i] += (pGrid[i] - mGrid[i]) * slew : pGrid[i];
			_grid[i] = mGrid[i];
		}
		
		int _slices = bSlices;

		mMirrorDistance = (morphing) ? mMirrorDistance += (mirrorDistance - mMirrorDistance) * slew : mirrorDistance;
		
		
		
		
		// - - - - - - - - - - - - - - - - - - - - - - -
		
		int _i = 0;			// current item
		boolean goon = true;	// go on
		
		float ratio = (float) texture_width / (float) texture_height;
		float _x = (float) Math.floor(ratio / _grid[0]) * -_grid[0];
		float _y = (float) Math.floor(1.f / _grid[1]) * -_grid[1];
		float borderx = -_x;
		float bordery = -_y;
		int row = 0;
		float flip = 1.f;
		
		while(goon) {

			drawPaisley(_i, _pos[0] + _x,_pos[1] + _y, _size[0], _size[1]*flip, _slices);
			
			
			// adapt the x and y values to the paisley grid
			// every second row the elements are positioned in between the top row grid
			_x+= _grid[0]*2;
			if(_x > borderx) {
				row++;
				_x = (row%2 == 0) ? -borderx : -borderx + _grid[0];
				flip = (row%2 ==0) ? 1f : -1f;
				_y += _grid[1];
				if(_y > bordery) goon = false;
			}

			_i++;
			if(_i >= pNumMax) goon = false;
		}

		
//		if(debug) post("drawimmediate");
		// call drawimmediate, to execute drawing of sketch object
		
//		try {
			sketch.call("drawimmediate");	
//		} catch(Exception e) {
//			if(debug) post("drawimmediate error: "+e);
//		}
		
		
	}


	
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === === paisley functions ==== === === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	/* create the master pattern */
	private void createMasterPaisley() {
		
		// bezier 1
//		masterPattern[0][0] = new float[] { 0.f, 0.f, 0.f };
//		masterPattern[0][1] = new float[] { 1.f, 0.f, 0.f };
//		masterPattern[0][2] = new float[] { -1.f, 1.f, 0.f };
//		masterPattern[0][3] = new float[] { 0.f, 1.f, 0.f };
		
		eNum = 4;
		
		int e = -1;
		
		e++; 
		// leaf
		eAnchor[e] = 1;
		masterWidth[e] = 0.1f;
		masterpattern(e,4,-0.5f, -0.43f, -0.3f, 0.08f, 0.f, -0.15f, 0.14f, 0.15f, 0, 0, 0, 0);
		
		e++;
		// leaf 2
		eAnchor[e] = 0;
		masterWidth[e] = 0.07f;
		masterpattern(e,4,-0.51f, -0.13f, -0.3f, 0.67f, 0.08f, 0.22f, -0.1f, 0.16f, 0, 0, 0, 0);
		
		e++; 
		// bottom curve
		eAnchor[e] = 2;
		masterWidth[e] = 0.05f;
		masterpattern(e,4,0.f,0.f, 0.f, 1.24f, -0.2f, 0.48f, -0.4f, 0.73f, 0, 0, 0, 0);
		
		
		e++;
		// drop
		eAnchor[e] = 0;
		masterWidth[e] = 0.07f;
		masterpattern(e,3,-0.06f, -0.24f, 0.27f, -0.05f, 0.29f, -0.18f, 0, 0, 0, 0, 0, 0);
		
		
		e++;
		// drop
//		eAnchor[e++] = 3;
//		masterWidth[e] = 0.07f;
//		masterpattern(e,3,-0.06f, 1.36f, -0.27f, 1.55f, 0.29f, 1.42f, 0, 0, 0, 0, 0, 0);
	}
	
	
	/* randomize bezier curves and widths */
	public void randomize(float f, float r) {
		randomizeElements(f,r);
		randomizeWidth(f,r);
	}
	
	
	/* initialize pattern bezier elements with random values */
	public void randomizeElements(float f, float r) {
		
		for(int p=0; p<pNumMax; p++) {
			for(int e=0; e<eNumMax; e++) {
				
				// randomize bezier point noise
				for(int b=0; b<bNum; b++) {
					// set noise to 0 if element is not active (<eNum)
					if(dice(f)) patternNoise[p][e][b][0] = (e<eNum) ? (float) Math.random()*2 - 1.f : 0.f;
					else if(dice(r)) patternNoise[p][e][b][0] = 0.f;
					if(dice(f)) patternNoise[p][e][b][1] = (e<eNum) ? (float) Math.random()*2 - 1.f : 0.f;
					else if(dice(r)) patternNoise[p][e][b][1] = 0.f;
					if(dice(f)) patternNoise[p][e][b][2] = 0.f;
					else if(dice(r)) patternNoise[p][e][b][2] = 0.f;
				}
			}
		}

	}
	
	/* chance function. return true if random value within threshold */
	private boolean dice(double f) {
		
		double v = Math.random();
		if(v < f) return true;
		else return false;
	}



	public void randomizeWidth(float f, float r) {
		for(int p=0; p<pNumMax; p++) {
			for(int e=0; e<eNumMax; e++) {
				// randomize width noise
				if(dice(f)) widthNoise[p][e] = (float) Math.random()*2 - 1.f;
				else if(dice(r)) widthNoise[p][e] = 0.f;
			}
		}
	}

	
	public void drawPaisley(int _i, float _x, float _y, float _sx, float _sy, int _slices) {
		
		// mWidthNoise
		//mPatternNoise
		
		// morphing
		float _pNoise[][][] = new float[eNum][bNum][3];
		float _wNoise[] = new float[eNum];

		for(int e=0; e<eNum; e++) {
			for(int b=0; b<bNum; b++) {
				for (int i = 0; i < 2; i++) {		// only need x and y values
					mPatternNoise[_i][e][b][i] = (morphing) ? mPatternNoise[_i][e][b][i] += (patternNoise[_i][e][b][i] - mPatternNoise[_i][e][b][i]) * slew : patternNoise[_i][e][b][i];
					_pNoise[e][b][i] = mPatternNoise[_i][e][b][i];
				}
			}
			mWidthNoise[_i][e] = (morphing) ? mWidthNoise[_i][e] += (widthNoise[_i][e] - mWidthNoise[_i][e]) * slew : widthNoise[_i][e];
			_wNoise[e] = mWidthNoise[_i][e];
		}
		
		float[] firstpoint = new float[] {0,0,0};
		float[] lastpoint = new float[] {0,0,0};
		
		for(int e=0; e<eNum; e++) {
			// add noise to masterpattern bezier curve
			float[][] masterDeviation = sumVector(masterPattern[e], scaleVector(_pNoise[e],bNoiseWeight), ePoints[e]);
			
			
			switch(eAnchor[e]) {
			case 1: 	// first point is fixed, revert to masterpattern x/y
				masterDeviation[0][0] = masterPattern[e][0][0];
				masterDeviation[0][1] = masterPattern[e][0][1];
				break;
			case 2:		// first point is attached to startpoint of last element, all other points are relative
				masterDeviation[0] = firstpoint;
				for(int i=1; i<ePoints[e]; i++) {
					masterDeviation[i][0]+=firstpoint[0];
					masterDeviation[i][1]+=firstpoint[1];
				}
				break;
			case 3: 	// first point is attached to endpoint of last element, all other points are relative
				masterDeviation[0] = lastpoint;
				for(int i=1; i<ePoints[e]; i++) {
					masterDeviation[i][0]+=lastpoint[0];
					masterDeviation[i][1]+=lastpoint[1];
				}
				break;
			}
			
			firstpoint = masterDeviation[0];
			lastpoint = (ePoints[e]==3) ? masterDeviation[0] : masterDeviation[ePoints[e]-1];
			
			// draw bezier curve
			drawBezier(_x, _y, masterDeviation, masterWidth[e] + _wNoise[e]*wNoiseWeight, _slices, _sx, _sy);
		}
	}
	
	
	/* draw n-grade bezier with strokewidth that converges at endpoints */
	public void drawBezier(float x, float y, float[][] points, float w, int slices, float _sx, float _sy) {
		
		// TODO calculate slices based on curve length ?
		float _add = 1.f / (float) (slices);			// defines segments of curve
		
		/* calculate all bezier points and their offset vectors beforehand */
		
		float[][] bezierPoint = new float[slices+1][3];
		float[][] offsetVector = new float[slices+1][3];
		
		// calculate all points along the bezier curve
		for(int _t=0; _t<=slices; _t++) {
			float t = _t*_add;
			bezierPoint[_t] = P(t, points);	
		}
		
		// calculate vector perpendicular to the line connecting the previous and next bezier point
		for(int _t=0; _t<=slices; _t++) {
			if(_t==0) offsetVector[_t] = perpendicularVector(bezierPoint[_t],bezierPoint[_t+1]);
			else if(_t==slices) offsetVector[_t] = perpendicularVector(bezierPoint[_t-1],bezierPoint[_t]);
			else offsetVector[_t] = perpendicularVector(bezierPoint[_t-1],bezierPoint[_t+1]);
		}
		
		// calculate offsetPoints of bands before drawing.
		float _scw = w;
		float[][][] bezierOffset = new float[slices+1][2][3];	// hold offset points on both sides of curve
		
		
		
		
		float _mirror = mMirrorDistance * _sx;
		float _dir = 1;							// mirror x values in opposite direction
		
		// mirror the bezier
		for(int m=0; m<2; m++) {
			if(m==1) _dir *= -1;

			for(int _t=0; _t<=slices; _t++) {
				float t = _t*_add;
				_scw = restrict(halfCircle(t) * w,0.01f,100f);
				float[] _offOutside = scaleVector(offsetVector[_t], _scw);
				float[] _offInside = scaleVector(offsetVector[_t], -_scw);

				bezierOffset[_t][0][0] = x + _mirror*_dir + bezierPoint[_t][0]*_sx*_dir + _offOutside[0]*_sx*_dir;
				bezierOffset[_t][1][0] = x + _mirror*_dir + bezierPoint[_t][0]*_sx*_dir + _offInside[0]*_sx*_dir;
				bezierOffset[_t][0][1] = y + bezierPoint[_t][1]*_sy + _offOutside[1]*_sy;
				bezierOffset[_t][1][1] = y + bezierPoint[_t][1]*_sy + _offInside[1]*_sy;
				bezierOffset[_t][0][2] = 0;
				bezierOffset[_t][1][2] = 0;
			}

			// --- - - - - - - - - --- - -draw BEZIER - - - -- - - -- - -- - -- -- -- - ---- 

			sketch.call("glcolor", pColor());

			sketch.call("glbegin", "tri_strip");

			for(int i=0; i<=slices; i++) {
				sketch.call("glvertex", new Atom[]{	Atom.newAtom(bezierOffset[i][0][0]),
						Atom.newAtom(bezierOffset[i][0][1]),
						Atom.newAtom(bezierOffset[i][0][2])});
				// if(i!=0 && i!=slices)  // don't draw the first and last point twice!
				sketch.call("glvertex", new Atom[]{	Atom.newAtom(bezierOffset[i][1][0]),
						Atom.newAtom(bezierOffset[i][1][1]),
						Atom.newAtom(bezierOffset[i][1][2])});
			}

			sketch.call("glend");

			// --- - - - - - - - - --- - - --end drawing -- - - -- - -- - -- -- -- - ---- -
		}
		
	}
	



	/* return color in atom format */
	private Atom[] pColor() {
		return new Atom[]{Atom.newAtom(1.0f),Atom.newAtom(1.0f),Atom.newAtom(1.0f)};
	}
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === ===  bezier functions  === === === === === === === 
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	/* Bezier : Computes factorial */
	// bezier computation based on
	// http://html5tutorial.com/how-to-draw-n-grade-bezier-curve-with-canvas-api/
	private float fact(float k) {
		if (k == 0 || k == 1) {
			return 1;
		} else {
			return k * fact(k - 1);
		}
	}

	/*
	 * Bezier : Computes Bernstain
	 * 
	 * @param {Integer} i - the i-th index
	 * 
	 * @param {Integer} n - the total number of points
	 * 
	 * @param {Number} t - the value of parameter t , between 0 and 1
	 */
	private float B(int i, int n, float t) {
		return (float) (fact(n) / (fact(i) * fact(n - i)) * Math.pow(t, i) * Math
				.pow(1 - t, n - i));
	}

	/*
	 * Computes a point's coordinates for a value of t
	 * 
	 * @param {Number} t - a value between o and 1
	 * 
	 * @param {Array} points - an {Array} of [x,y] coordinates. The initial
	 * points
	 */
	private float[] P(float t, float[][] points) {
		float[] r = { 0, 0, 0 };
		int n = points.length - 1;
		for (int i = 0; i <= n; i++) {
			r[0] += points[i][0] * B(i, n, t);
			r[1] += points[i][1] * B(i, n, t);
			r[2] += points[i][2] * B(i, n, t);
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
	
	private float[][] scaleVector(float[][] v, float s) {
		int a1 = v.length;
		int a2 = v[0].length;
		float[][] newVector = new float[a1][a2];
		for(int i=0; i<a1; i++) {
			for(int j=0; j<a2; j++) {
				newVector[i][j] = v[i][j] * s;
			}
		}
		return newVector;
	}
	
	private float[][] sumVector(float[][] v1, float[][] v2) {
		int a1 = v1.length;
		return sumVector(v1, v2, a1);
	}
	
	/* no indicates the amount of array positions the sum vector should have */
	private float[][] sumVector(float[][] v1, float[][] v2, int no) {
		int a1 = no;
		if(no==3) a1 = 4;
		int a2 = v1[0].length;
		float[][] sumVector = new float[a1][a2];
		for(int i=0; i<no; i++) {
			for(int j=0; j<a2; j++) {
				sumVector[i][j] = v1[i][j] + v2[i][j];
			}
		}
		if(no==3) {
			for(int j=0; j<a2; j++) sumVector[3][j] = v1[0][j] + v2[0][j];
		}
		return sumVector;
	}
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === == helpful little functions == === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	private float restrict(float v) {
		return restrict(v, -1.f, 1.f);
	}

	private float restrict(float v, float min, float max) {
		if (v < min)
			return min;
		if (v > max)
			return max;
		return v;
	}
	
	
	
	/* translate from 0 to 1, to a half circle
	 * 0.0 = 0
	 * 0.5 = 1
	 * 1.0 = 0
	 */
	private float halfCircle(float t) {
		return (float) Math.sin(t * 3.1415927f);
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
		lineSmooth = (v == 1) ? true : false;
	}
	
	
	/* === === === === === === === paisley === === === === === === === === */

	/* set center position of rainbow */
	public void position(float x, float y, float z) {
		pPos[0] = x;
		pPos[1] = y;
		pPos[2] = z;
	}
	
	/* set main scale of rainbow */
	public void size(float x, float y, float z) {
		pSize[0] = x;
		pSize[1] = y;
		pSize[2] = z;
	}
	
	public void grid(float x, float y) {
		pGrid[0] = x;
		pGrid[1] = y;
	}
	
	public void debug(int v) {
		debug = (v==1) ? true : false;
	}
	
	
	/* toggle morphing */
	public void morph(int v) {
		morphing = (v == 1) ? true : false;
	}

	/* set morphing speed */
	public void morphspeed(float v) {
		slew = (v > 0) ? v : 0.01f;
	}
	
	public void animate(int v) {
		animation = (v==1) ? true : false;
	}
	
	public void noisefactor(float v) {
		noiseFactor = v;
	}
	
	public void returnfactor(float v) {
		returnFactor = v;
	}
	
	/* set the weight of the pattern noise */
	public void noiseweight(float v1, float v2) {
		bNoiseWeight = v1;
		wNoiseWeight = v2;
	}
	
	public void masterwidth(float v) {
		for(int e=0; e<eNumMax; e++ ) {
			masterWidth[e] = v;
		}
	}
	
	/* resolution of bezier curve */
	public void slices(int v) {
		bSlices = (v>2) ? v : 2;
	}
	
	public void mirrordistance(float v) {
		mirrorDistance = v;
	}
	
	/* input bezier curve points */
	public void masterpattern(int e, int no, float x1, float y1, float x2, float y2, float x3, float y3,
			float x4, float y4, float x5, float y5, float x6, float y6) {
		ePoints[e] = no;
		masterPattern[e][0] = new float[] { x1, y1, 0.f };
		masterPattern[e][1] = new float[] { x2, y2, 0.f };
		masterPattern[e][2] = new float[] { x3, y3, 0.f };
		masterPattern[e][3] = new float[] { x4, y4, 0.f };
		masterPattern[e][4] = new float[] { x5, y5, 0.f };
		masterPattern[e][5] = new float[] { x6, y6, 0.f };
	}
	
	
	/* === === === === === === === GUI in max patch === === === === === === === === */
	

	public void outputVariables() {
		outputAllVariables();
	}
	
	
	private void outputAllVariables() {

		// display parameters
		outlet(1,"script send gui_depthenable set "+sketchDepthEnable);
		outlet(1,"script send gui_antialias set "+sketchAntialias);
		outlet(1,"script send gui_fsaa set "+sketchFsaa);
		outlet(1,"script send gui_linesmooth set "+ (lineSmooth ? 1 : 0) );
		outlet(1,"script send gui_blendenable set "+sketchBlendEnable);
		outlet(1,"script send gui_blendmode1 set "+sketchBlendMode1);
		outlet(1,"script send gui_blendmode2 set "+sketchBlendMode2);
		
		outlet(1,"script send gui_animate set "+(animation ? 1 : 0));
		outlet(1,"script send gui_noisefactor set "+noiseFactor);
		outlet(1,"script send gui_returnfactor set "+returnFactor);
		outlet(1,"script send gui_morphspeed set "+slew);
		outlet(1,"script send gui_morph set "+ (morphing ? 1 : 0));

		outlet(1,"script send gui_size0 set "+pSize[0]);
		outlet(1,"script send gui_sizex set "+pSize[0]);
		outlet(1,"script send gui_sizey set "+pSize[1]);
		outlet(1,"script send gui_sizez set "+pSize[2]);
		outlet(1, "script send gui_size set size " + pSize[0] + " " + pSize[1] + " " + pSize[2]);
		outlet(1,"script send gui_positionx set "+pPos[0]);
		outlet(1,"script send gui_positiony set "+pPos[1]);
		outlet(1,"script send gui_positionz set "+pPos[2]);
		outlet(1, "script send gui_position set position " + pPos[0] + " " + pPos[1] + " " + pPos[2]);
		
		outlet(1,"script send gui_gridx set "+pGrid[0]);
		outlet(1,"script send gui_gridy set "+pGrid[1]);
		outlet(1, "script send gui_grid set grid " + pGrid[0] + " " + pGrid[1]);
		
		outlet(1,"script send gui_noiseweight1 set "+bNoiseWeight);
		outlet(1,"script send gui_noiseweight2 set "+wNoiseWeight);
		outlet(1,"script send gui_noiseweight set noiseweight "+bNoiseWeight + " " +wNoiseWeight);
		
		outlet(1,"script send gui_slices set "+bSlices);
		outlet(1,"script send gui_masterwidth set "+masterWidth[0]);
		outlet(1,"script send gui_mirrordistance set "+mirrorDistance);
		
		
	}
	

	
	
	
	
}
