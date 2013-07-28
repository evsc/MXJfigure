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
	private int texture_width = 1600; 		
	private int texture_height = 1200;
	
	
	private float[] bgColor = { 0.f, 0.f, 0.f, 1.f };
	private float[] bgColorFade1 = { 0,0,0,1 };
	private float[] bgColorFade2 = { 1,1,1,1 };
	private boolean bgFade = false;
	private boolean bgFading = false;		// set true while fading
	private int bgFadeTime = 25;
	private int bgFadeCounter = 0;
	
	// grid
	private float[] pGrid = { 0.6f, 0.6f };
	private float[] pPos = {0,0,0};			// center position
	private float[] pSize = {1,1,1};		// size of all (grid + pattern size)
	private float[] pScale = {2,2,2};		// size of the pattern
	private float mirrorDistance = 0.0f;
	private int displayMode = 0; 			// grid ... 0, twin ... 1, one ... 2
	private boolean flipTwin = false;		// flip pattern of twin mode
	private float cutOff = 0.f;				// lower part of grid that doesn't get displayed
	
	// 
	private int pNumMax = 96;				// maximum number of paisley patterns
	private int rows = 6;
	private int cols = 16;
	
	// pattern parameter
	private float[] mainColor = { 1.f, 1.f, 1.f, 1.f };
	private float[][] patternColor;	
	private float[][][] masterPattern;		// the master pattern template, 
											// every displayed pattern will be a slight variation of this one
	
	float[][][][] masterBezierOffset;		// hold calculated bezier points for easy drawing of master pattern
	
	private int eNum = 4;					// number of elements that constitute one paisley pattern
	private int eNumMax = 5;
	private boolean _mirrored = true;
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
	private boolean precalc = true;
	private int[] theActive = { 39, 56, 58, 70, 72 };
	private boolean animationQueue = false;	// only animate one after the other
	private int animationQueuePointer = 0;	// which of the 5 patterns
	private int animationQueueCounter = 0;
	private int animationQueueTime = 20; 
	private boolean[] activated;			// which pattern is allowed to move
	private boolean[] displayed;			// which pattern are displayed
	private float noiseFactor = 1.f;		// how likely a noisevalue will be changed
	private float returnFactor = 10.f;		// how likely the noisevalue will be set to 0
	private float slew = 0.1f; 				// morphing speed (position, size, grid, ...)
	private float animationSpeed = 0.1f;	// morphing speed of pattern morphing
	private float[] mPos = { 0, 0, 0 };
	private float[] mSize = { 1, 1, 2 };
	private float[] mScale = { 2, 2, 2 };
	private float[] mGrid = { 0.4f, 0.4f };
	private float mMirrorDistance = 0.1f;
	private float scaleMult = 0.1f;			// scale down size, as 100% would cover whole screen
	private float[][][][] mPatternNoise;	// morphing the noise deviation
	private float[][] mWidthNoise;			// morphing noise width deviation
	private float[] mMasterWidth;			
	private float mBNoiseWeight = 0.05f;
	private float mWNoiseWeight = 0.1f;
	
	
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
		mMasterWidth = new float[eNumMax];
		widthNoise = new float[pNumMax][eNumMax];
		mWidthNoise = new float[pNumMax][eNumMax];
		ePoints = new int[eNumMax];
		eAnchor = new int[eNumMax];
		activated = new boolean[pNumMax];
		displayed = new boolean[pNumMax];
		patternColor = new float[pNumMax][3];
		
		for(int e=0; e<eNumMax; e++ ) {
			masterWidth[e] = 0.05f;
			mMasterWidth[e] = 0.05f;
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
			activated[p] = false;
			displayed[p] = true;
			for(int c=0; c<3; c++) patternColor[p][c] = 1.f;
		}
		
		for(int i=0; i<5; i++) activate(theActive[i],1);
		
		createMasterPaisley();
		randomize(0.1f,0.f);
		display(displayMode);
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
		
		if(bgFading) {
			bgFadeCounter++;
			float fadeStep = bgFadeCounter * 1.f/(float) bgFadeTime;
			if(debug) post("fadestep "+fadeStep+"  counter "+bgFadeCounter);
			bgColor[0] = bgColorFade1[0] + (bgColorFade2[0]-bgColorFade1[0]) * fadeStep;
			bgColor[1] = bgColorFade1[1] + (bgColorFade2[1]-bgColorFade1[1]) * fadeStep;
			bgColor[2] = bgColorFade1[2] + (bgColorFade2[2]-bgColorFade1[2]) * fadeStep;
			
			if(bgFadeCounter >= bgFadeTime) {
				// reached end of fade
				bgColor[0] = bgColorFade2[0];
				bgColor[1] = bgColorFade2[1];
				bgColor[2] = bgColorFade2[2];
				bgFading = false;
			}
		}

//		if(debug) post("reset");
		sketch.call("reset");
		
		sketch.call("glclearcolor",
				new Atom[] { 	Atom.newAtom(bgColor[0]), 
								Atom.newAtom(bgColor[1]), 
								Atom.newAtom(bgColor[2]), 
								Atom.newAtom(bgColor[3]) });
		sketch.call("glclear");
		
//		if(debug) post("line_smooth");
		if (lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");
		
		// create local variables, to avoid conflict when life-updating
		// variables while rendering
		
		float _size[] = new float[2];		// size from origin
		float _scale[] = new float[2];
		float _pos[] = new float[2];		// origin position
		float _grid[] = new float[2];
		
		for (int i = 0; i < 2; i++) {		// only need x and y values
			mSize[i] = (morphing) ? mSize[i] += (pSize[i] - mSize[i]) * slew : pSize[i];
			_size[i] = mSize[i];
			mScale[i] = (morphing) ? mScale[i] += (pScale[i] - mScale[i]) * slew : pScale[i];
			_scale[i] = mScale[i] * scaleMult;
			mPos[i] = (morphing) ? mPos[i] += (pPos[i] - mPos[i]) * slew : pPos[i];
			_pos[i] = mPos[i];
			mGrid[i] = (morphing) ? mGrid[i] += (pGrid[i] - mGrid[i]) * slew : pGrid[i];
			_grid[i] = mGrid[i];
		}
		
		int _slices = bSlices;

		mMirrorDistance = (morphing) ? mMirrorDistance += (mirrorDistance - mMirrorDistance) * slew : mirrorDistance;
		
		
		
	
		
		
		// new drawing loop, create 16x9 grid and only draw if item is within frame
		float borderx = -_grid[0] * cols/2.f;
		float bordery = _grid[1] * rows/2.f;
		float wratio = 0.1f + (float) texture_width / (float) texture_height;
		float hratio = 0.1f + 1.f;		// add 0.1f to avoid on/off flickr on edges
		float flip = 1.f;
		int count=0;					// to calculate how many patterns are actually displayed
		for(int r=0; r<rows; r++) {
			for(int c=0; c<cols; c++) {
				int id= r*cols +c;
				if(displayed[id]) {
					float _x = borderx + c * _grid[0];
					_x += (r%2 == 0) ? _grid[0]/2.f : 0;	// each second is shifted to the right
					float _y = bordery - r * _grid[1];
					flip = (r%2 ==0) ? 1f : -1f;
					_x *= _size[0];
					_y *= _size[1];
					if(_x >= -wratio && _x <= wratio && _y >= -hratio+cutOff*2 && _y <= hratio) {
						count++;
						if(displayMode==1) {
							// twin
							flip *= flipTwin ? -1.f : 1.f;
							drawPaisley(id, (_pos[0] + _x),(_pos[1] + _y-_grid[1]*_size[1]), _scale[0]*_size[0], _scale[1]*flip*_size[1], activated[id] && precalc);	
							drawPaisley(id, (_pos[0] + _x),(_pos[1] + _y), _scale[0]*_size[0], _scale[1]*flip*-1*_size[1], activated[id] && precalc);	
						} else {
							drawPaisley(id, (_pos[0] + _x),(_pos[1] + _y), _scale[0]*_size[0], _scale[1]*flip*_size[1], activated[id] && precalc);	
						}
					}
				}
			}
		}
		
//		post("drawing "+count+" patterns");
		
//		if(debug) post("drawimmediate");
		// call drawimmediate, to execute drawing of sketch object
		
//		try {
			sketch.call("drawimmediate");	
//		} catch(Exception e) {
//			if(debug) post("drawimmediate error: "+e);
//		}
//		if(debug) post("after drawimmediate");
		
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
//		masterpattern(e,4,-0.5f, -0.63f, -0.3f, -0.12f, 0.f, -0.35f, 0.14f, -0.05f, 0, 0, 0, 0);
		
		e++;
		// leaf 2
		eAnchor[e] = 0;
		masterWidth[e] = 0.07f;
		masterpattern(e,4,-0.51f, -0.13f, -0.3f, 0.67f, 0.08f, 0.22f, -0.1f, 0.16f, 0, 0, 0, 0);
//		masterpattern(e,4,-0.51f, -0.33f, -0.3f, 0.47f, 0.08f, 0.02f, -0.1f, -0.04f, 0, 0, 0, 0);
		
		e++; 
		// bottom curve
		eAnchor[e] = 2;
		masterWidth[e] = 0.05f;
		masterpattern(e,4,0.f,0.f, 0.f, 1.24f, -0.2f, 0.48f, -0.4f, 0.73f, 0, 0, 0, 0);
//		masterpattern(e,4,0.f,-0.2f, 0.f, 1.04f, -0.2f, 0.28f, -0.4f, 0.53f, 0, 0, 0, 0);
		
		
		e++;
		// drop
		eAnchor[e] = 0;
		masterWidth[e] = 0.07f;
		masterpattern(e,3,-0.06f, -0.24f, 0.27f, -0.05f, 0.29f, -0.18f, 0, 0, 0, 0, 0, 0);
//		masterpattern(e,3,-0.06f, -0.44f, 0.27f, -0.25f, 0.29f, -0.38f, 0, 0, 0, 0, 0, 0);
		
		/* special anchoring of element 
			*  0 = no effect
			*  1 = keep first point fixed
			*  2 = first point is relative to startpoint of previous element
			*  3 = first point is relative to endpoint of previous element
			*/
		
		//y shift
		// masterPattern = new float[eNumMax][bNum][3]
		for(int en=0; en<eNumMax; en++) {
			if(en==0 || en==3 || en==1) {
				for(int b=0; b<ePoints[en]; b++) {
					if(en==2 && b==0) {
						
					} else {
						masterPattern[en][b][1] = masterPattern[en][b][1] - 0.1f;
					}
					masterPattern[en][b][0] += 0.5f;
				}
			}
		}
		
		
		calculateMasterPattern();	
		
		
	}
	
	// calculate bezier curves for master
	public void calculateMasterPattern() {
		
		// global variable with resulting bezier points
		masterBezierOffset = new float[eNumMax][bSlices+1][2][3];
		
		float[] firstpoint = new float[] {0,0,0};
		float[] lastpoint = new float[] {0,0,0};
		
		float _add = 1.f / (float) (bSlices);			// defines segments of curve
		
		for(int e=0; e<eNum; e++) {
			
			float _empty[][] = new float[bNum][3];
			float[][] masterDeviation = sumVector(masterPattern[e], _empty, ePoints[e]);
			
			
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


			/* calculate all bezier points and their offset vectors beforehand */
			float[][] bezierPoint = new float[bSlices+1][3];
			float[][] offsetVector = new float[bSlices+1][3];

			// calculate all points along the bezier curve
			for(int _t=0; _t<=bSlices; _t++) {
				float t = _t*_add;
				bezierPoint[_t] = P(t, masterDeviation);	
			}

			// calculate vector perpendicular to the line connecting the previous and next bezier point
			for(int _t=0; _t<=bSlices; _t++) {
				if(_t==0) offsetVector[_t] = perpendicularVector(bezierPoint[_t],bezierPoint[_t+1]);
				else if(_t==bSlices) offsetVector[_t] = perpendicularVector(bezierPoint[_t-1],bezierPoint[_t]);
				else offsetVector[_t] = perpendicularVector(bezierPoint[_t-1],bezierPoint[_t+1]);
			}

			// calculate offsetPoints before drawing.
			float _scw = masterWidth[e];

			for(int _t=0; _t<=bSlices; _t++) {
				float t = _t*_add;
				_scw = restrict(halfCircle(t) * masterWidth[e],0.01f,100f);
				float[] _offOutside = scaleVector(offsetVector[_t], _scw);
				float[] _offInside = scaleVector(offsetVector[_t], -_scw);

				masterBezierOffset[e][_t][0][0] = bezierPoint[_t][0] + _offOutside[0];
				masterBezierOffset[e][_t][1][0] = bezierPoint[_t][0] + _offInside[0];
				masterBezierOffset[e][_t][0][1] = bezierPoint[_t][1] + _offOutside[1];
				masterBezierOffset[e][_t][1][1] = bezierPoint[_t][1] + _offInside[1];
				masterBezierOffset[e][_t][0][2] = -0.5f;
				masterBezierOffset[e][_t][1][2] = -0.5f;
			}
		}
	}
	
	
	/* randomize bezier curves and widths */
	public void randomize(float f, float r) {
		if(animationQueue) {
			animationQueueCounter++;
			if(animationQueueCounter > animationQueueTime) {
				animationQueueCounter = 0;
				animationQueuePointer++;
				if(animationQueuePointer>4) {
					animationQueuePointer = 0;
				}
			}
		}
		randomizeElements(f,r);
		randomizeWidth(f,r);
	}
	
	
	public void clearNoise() {
		randomizeElements(0,1,true);
		randomizeWidth(0,1,true);
	}
	
	public void randomizeElements(float f, float r) {
		randomizeElements(f,r,false);
	}
	
	/* initialize pattern bezier elements with random values */
	public void randomizeElements(float f, float r, boolean all) {
		
		
		for(int p=0; p<pNumMax; p++) {
			if(activated[p] || all) {
				if(all || !animationQueue || (animationQueue && p==theActive[animationQueuePointer])) {
//					post("p "+p+"   animationQueue "+animationQueue+"  animationQueuePointer "+animationQueuePointer+
//							"   theActive[animationQueuePointer]) "+theActive[animationQueuePointer]);			
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
		}

	}
	
	/* chance function. return true if random value within threshold */
	private boolean dice(double f) {
		
		double v = Math.random();
		if(v < f) return true;
		else return false;
	}

	public void randomizeWidth(float f, float r) {
		randomizeWidth(f,r,false);
	}

	public void randomizeWidth(float f, float r, boolean all) {
		for(int p=0; p<pNumMax; p++) {
			if(activated[p] || all) {
				if(all || !animationQueue || (animationQueue && p==theActive[animationQueuePointer])) {
					for(int e=0; e<eNumMax; e++) {
						// randomize width noise
						if(dice(f)) widthNoise[p][e] = (float) Math.random()*2 - 1.f;
						else if(dice(r)) widthNoise[p][e] = 0.f;
					}
				}
			}
		}
	}
	
	
	// function that draws already calculated pattern
	public void drawPaisley(int _i, float _x, float _y, float _sx, float _sy, boolean calc){
		
		if(calc) renderDrawPaisley(_i, _x, _y, _sx, _sy);
		else {
			float _mirror = mMirrorDistance * _sx;

			for(int e=0; e<eNum; e++) {

				float _dir = 1;							// mirror x values in opposite direction
				// mirror the bezier
				for(int m=0; m<2; m++) {
					if(m==0 || (m==1 && _mirrored)) {
						if(m==1) _dir *= -1;

						// --- - - - - - - - - --- - -draw BEZIER - - - -- - - -- - -- - -- -- -- - ---- 
						sketch.call("glcolor", pColor());		// same color for all patterns

						sketch.call("glbegin", "tri_strip");

						for(int i=0; i<=bSlices; i++) {
							// x + _mirror*_dir + bezierPoint[_t][0]*_sx*_dir + _offOutside[0]*_sx*_dir;
							// y + bezierPoint[_t][1]*_sy + _offOutside[1]*_sy
							float x1 = _x + _mirror*_dir + masterBezierOffset[e][i][0][0]*_sx*_dir;
							float y1 = _y + masterBezierOffset[e][i][0][1]*_sy;
							float z1 = 0;
							float x2 = _x + _mirror*_dir + masterBezierOffset[e][i][1][0]*_sx*_dir;
							float y2 = _y + masterBezierOffset[e][i][1][1]*_sy;
							float z2 = 0;

							sketch.call("glvertex", new Atom[]{	Atom.newAtom(x1), Atom.newAtom(y1), Atom.newAtom(z1)});
							sketch.call("glvertex", new Atom[]{	Atom.newAtom(x2), Atom.newAtom(y2), Atom.newAtom(z2)});
						}

						sketch.call("glend");
					}
					// --- - - - - - - - - --- - - --end drawing -- - - -- - -- - -- -- -- - ---- -
				}

			}
		}
	}
	
	
	// function that needs to renders bezier points
	public void renderDrawPaisley(int _i, float _x, float _y, float _sx, float _sy) {
		
		// mWidthNoise
		//mPatternNoise
		
		// morphing
		float _pNoise[][][] = new float[eNum][bNum][3];
		float _wNoise[] = new float[eNum];

		for(int e=0; e<eNum; e++) {
			for(int b=0; b<bNum; b++) {
				for (int i = 0; i < 2; i++) {		// only need x and y values
					if(animation) mPatternNoise[_i][e][b][i] += (patternNoise[_i][e][b][i] - mPatternNoise[_i][e][b][i]) * animationSpeed;
					_pNoise[e][b][i] = mPatternNoise[_i][e][b][i];
				}
			}
			if(animation) mWidthNoise[_i][e] += (widthNoise[_i][e] - mWidthNoise[_i][e]) * animationSpeed;
			_wNoise[e] = mWidthNoise[_i][e];
		}
		
		
		mWNoiseWeight = (morphing) ? mWNoiseWeight += (wNoiseWeight-mWNoiseWeight) * slew : wNoiseWeight;
		mBNoiseWeight = (morphing) ? mBNoiseWeight += (bNoiseWeight-mBNoiseWeight) * slew : bNoiseWeight;
		
		float[] firstpoint = new float[] {0,0,0};
		float[] lastpoint = new float[] {0,0,0};
		
		for(int e=0; e<eNum; e++) {
			// add noise to masterpattern bezier curve
			float[][] masterDeviation = sumVector(masterPattern[e], scaleVector(_pNoise[e],mBNoiseWeight), ePoints[e]);
			
			mMasterWidth[e] = (morphing) ? mMasterWidth[e] += (masterWidth[e] - mMasterWidth[e]) * slew : masterWidth[e];
			float _mw = mMasterWidth[e];
		
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
			drawBezier(_i, _x, _y, masterDeviation, _mw + _wNoise[e]*mWNoiseWeight, _sx, _sy);
		}
	}
	
	
	/* draw n-grade bezier with strokewidth that converges at endpoints */
	public void drawBezier(int _i, float x, float y, float[][] points, float w, float _sx, float _sy) {
		
		int slices = bSlices;
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
		
		// calculate offsetPoints before drawing.
		float _scw = w;
		float[][][] bezierOffset = new float[slices+1][2][3];	// hold offset points on both sides of curve

		float _mirror = mMirrorDistance * _sx;
		float _dir = 1;							// mirror x values in opposite direction
		
		// mirror the bezier
		for(int m=0; m<2; m++) {
			if(m==0 || (m==1 && _mirrored)) {
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
					bezierOffset[_t][0][2] = -0.f;
					bezierOffset[_t][1][2] = -0.f;
				}
	
				// --- - - - - - - - - --- - -draw BEZIER - - - -- - - -- - -- - -- -- -- - ---- 
	
//				sketch.call("glcolor", pColor(_i));		// patterns individual color
				sketch.call("glcolor", pColor());		// same color for all patterns
	
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
			}
			// --- - - - - - - - - --- - - --end drawing -- - - -- - -- - -- -- -- - ---- -
		}
		
	}
	
	private Atom[] pColor() {
//		mainColor
		return new Atom[]{
				Atom.newAtom(mainColor[0]),
				Atom.newAtom(mainColor[1]),
				Atom.newAtom(mainColor[2])};
	}


	/* return color in atom format */
	private Atom[] pColor(int _i) {
		return new Atom[]{
				Atom.newAtom(patternColor[_i][0]),
				Atom.newAtom(patternColor[_i][1]),
				Atom.newAtom(patternColor[_i][2])};
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

	public void clearColor(float r, float g, float b) {
		if(debug) post("clearcolor() "+r+" "+g+" "+b);
		if(bgFade) {
			// set start color
			bgColorFade1[0] = bgColor[0];
			bgColorFade1[1] = bgColor[1];
			bgColorFade1[2] = bgColor[2];
			
			bgFading = true;		// turn on fading, so its processed during draw()
			bgFadeCounter = 0;
		}
		
		// set goal color
		bgColorFade2[0] = (float) r / 255.f;
		bgColorFade2[1] = (float) g / 255.f;
		bgColorFade2[2] = (float) b / 255.f;
		
		if(!bgFade) {
			bgColor[0] = bgColorFade2[0];
			bgColor[1] = bgColorFade2[1];
			bgColor[2] = bgColorFade2[2];
		}
	}

	public void fsaa(int v) {
		sketchFsaa = (v == 1) ? 1 : 0;
		sketch.setAttr("fsaa", sketchFsaa);
	}

	public void linesmooth(int v) {
		lineSmooth = (v == 1) ? true : false;
	}
	
	
	/* === === === === === === === paisley === === === === === === === === */
	
	public void bgfade(int v) {
		bgFade = (v==1) ? true : false;
	}
	
	public void bgfadetime(int v) {
		bgFadeTime = v;
	}

	/* set center position */
	public void position(float x, float y, float z) {
		pPos[0] = x;
		pPos[1] = y;
		pPos[2] = z;
	}
	
	/* set main scale */
	public void size(float x, float y, float z) {
		pSize[0] = x;
		pSize[1] = y;
		pSize[2] = z;
	}
	
	public void scale(float x, float y, float z) {
		pScale[0] = x;
		pScale[1] = y;
		pScale[2] = z;
	}
	
	public void grid(float x, float y) {
		pGrid[0] = x;
		pGrid[1] = y;
	}
	
	public void debug(int v) {
		debug = (v==1) ? true : false;
	}
	
	public void activate(int no, int v) {
		if(no >= 0 && no <= pNumMax) {
			activated[no] = (v==0) ? false : true;
		}
	}
	
	public void activeFive(int no1, int no2, int no3, int no4, int no5) {
		for(int p=0; p<pNumMax; p++) activated[p] = false;
		theActive = new int[] { no1, no2, no3, no4, no5 };
		for(int i=0; i<5; i++) activate(theActive[i],1);
	}
	
	public void display(int v) {
		displayMode = v;
		_mirrored = (v==2) ? false : true;
		if(displayMode==0) {
			for(int p=0; p<pNumMax; p++) displayed[p] = true;
		} else {
			for(int p=0; p<pNumMax; p++) displayed[p] = false;
			displayed[56] = true;
		}
	}
	
	/* toggle morphing */
	public void morph(int v) {
		morphing = (v == 1) ? true : false;
	}

	/* set morphing speed */
	public void morphspeed(float v) {
		slew = (v > 0) ? v : 0.001f;
	}
	
	public void animationspeed(float v) {
		animationSpeed = (v > 0) ? v : 0.001f;
	}
	
	public void animate(int v) {
		animation = (v==1) ? true : false;
	}
	
	public void calc(int v) {
		precalc = (v==1) ? true : false;
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
		calculateMasterPattern();	// recalculate master pattern bezier points
	}
	
	/* resolution of bezier curve */
	public void slices(int v) {
		bSlices = (v>2) ? v : 2;
		calculateMasterPattern();	// recalculate master pattern bezier points
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
	
	public void patterncolor(int _i, int r, int g, int b) {
		if(_i >= 0 && _i < pNumMax) {
			patternColor[_i][0] = (float) r / 255.f;
			patternColor[_i][1] = (float) g / 255.f;
			patternColor[_i][2] = (float) b / 255.f;
		}
	}
	
	public void maincolor(int r, int g, int b) {
		mainColor[0] = (float) r / 255.f; 
		mainColor[1] = (float) g / 255.f; 
		mainColor[2] = (float) b / 255.f; 
	}
	
	public void mirrored(int v) {
		_mirrored = (v==1) ? true : false;
	}
	
	public void animationqueue(int v) {
		animationQueue = (v==1) ? true : false;
	}
	
	public void queuetime(int v) {
		animationQueueTime = v;
	}
	
	public void cutoff(float v) {
		cutOff = v;
	}
	
	public void fliptwin(int v) {
		flipTwin = (v==1) ? true : false;
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
		outlet(1,"script send gui_animationspeed set "+animationSpeed);
		outlet(1,"script send gui_morph set "+ (morphing ? 1 : 0));

		outlet(1,"script send gui_size0 set "+pSize[0]);
//		outlet(1,"script send gui_sizex set "+pSize[0]);
//		outlet(1,"script send gui_sizey set "+pSize[1]);
//		outlet(1,"script send gui_sizez set "+pSize[2]);
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
		
		
		outlet(1,"script send gui_clearcolor set clearColor "
				+(int) (bgColor[0]*255)+" "+(int) (bgColor[1]*255)+" "+(int) (bgColor[2]*255));
		outlet(1,"script send gui_clearcolor1 bgcolor "+bgColor[0]+" "+bgColor[1]+
				" "+bgColor[2] + " 1.");
		outlet(1,"script send gui_bgfade set " + (bgFade ? 1 : 0));
		outlet(1,"script send gui_bgfadetime set "+ bgFadeTime);
		
		outlet(1,"script send gui_mirrored set "+(_mirrored ? 1 : 0));
		outlet(1,"script send gui_maincolor set maincolor "
				+(int) (mainColor[0]*255)+" "+(int) (mainColor[1]*255)+" "+(int) (mainColor[2]*255));
		outlet(1,"script send gui_maincolor1 bgcolor "+mainColor[0]+" "+mainColor[1]+
				" "+mainColor[2] + " 1.");
		
		outlet(1,"script send gui_queuetime set "+animationQueueTime);
		outlet(1,"script send gui_animationqueue set "+(animationQueue? 1: 0));
		
		outlet(1,"script send gui_active1 set "+theActive[0]);
		outlet(1,"script send gui_active2 set "+theActive[1]);
		outlet(1,"script send gui_active3 set "+theActive[2]);
		outlet(1,"script send gui_active4 set "+theActive[3]);
		outlet(1,"script send gui_active5 set "+theActive[4]);
		
		outlet(1,"script send gui_cutoff set "+cutOff);
		
	}
	

	
	//notifyDeleted is called by the Max application
	//when the user deletes your external from a Max patch
	//or closes a Max patch of which your Java extern
	//is a member.
	public void notifyDeleted()
	{
		// free max peers. otherwise these will persist for a while
		// until the garbage collector feels like cleaning up 
		texture.freePeer();
		sketch.freePeer();
	}
	
	
}
