/*
 * Figure de la terre - mxj objects
 * 
 * aurora borealis
 * main bezier curve
 * with multiple curtains across different angles
 * 
 */

import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;

public class FigureAuroraMulti extends MaxObject {

	// render context for jitter opengl objects
	private String context = "foo";
	private boolean debug = false;

	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	private int texture_width = 640; 	// TODO: include resize function
	private int texture_height = 480;

	// Bezier curve, aurora backbone
	private float[] aPos = { 0, 0, 0 };
	private float aBeta = 0.f; 			// birdview angle between outer bezier points
	private float[] aSize = { 1, 1, 1 };
	private int aNumMax = 3; 			// maximum number of aurora curtains
	private int aNum = 1; 				// number of active curtains
	private int abPointCount = 4; 		// main bezier point count, leave count
										// flexible just in case
	private float abPoint[][][]; 		// main bezier points
	private float abHeight[];		 	// z of inner bezier points defines height of
										// curve
	private float abAngle[]; 			// angle of curve, leaning right or left ?

	// aurora rays
	private int rayCount[];
	private float rayHeight[];
	private float rayWidth[];
	private int raySegments = 4; // color segmentation for rays (smaller =
									// better for cpu)
	private boolean rayColorByHeight = true;
	private float rayColorOffset[];
	private float rayColorMult[];
	private float[] aVanishingPoint = { 0, 100, 0 };

	// noise - main noise field
	private int noiseCount[];
	private int noiseF[]; // noise frequency
	private int noiseP[]; // noise pointer, for animating the noise
	private int noiseStep[];
	private boolean noiseAnimation = false;
	private int noiseEdge = 0;
	private boolean noiseSmoothEdge = false;
	private float[][][] noise; // curtains+1 / xyz / noiseCount
	private float noiseWeight[][]; // curtains+1 / xyz

	// morphing
	private boolean morphing = false;
	private float slew = 0.1f; // morphing speed
	private int mRayCount[];
	private float mPos[] = { 0, 0, 0 };
	private float mSize[] = { 1, 1, 1 };
	private float mbPoint[][][];
	private float[] mVanishingPoint = { 0, 100, 0 };
	// private float mbHeight = 0.5f;
	// private float mbAngle = 0.0f;
	private float mRayHeight[];

	// display parameters
	private boolean lineSmooth = false;
	private int sketchDepthEnable = 0;
	private int sketchAntialias = 0;
	private int sketchFsaa = 0;
	private int sketchBlendEnable = 0;

	/*
	 * === === === === === === === main functions === === === === === === ===
	 */

	/* instantiating mxj without argument */
	public FigureAuroraMulti() {
		bail("please provide render context as argument");
	}

	/* instantiate mxj with render context as argument */
	public FigureAuroraMulti(String rc) {

		context = rc; // render context
		declareIO(2, 1); // declare 2 inlets, 1 outlet of DataTypes.ALL

		// assist message for inlets and outlets (mouse hover)
		setInletAssist(new String[] { "bang to compute and draw",
				"input settings" });
		setOutletAssist(new String[] { "outputs jit.gl.texture object" });

		// instantiate Jitter sketch object
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", sketchDepthEnable);
		sketch.setAttr("blend_enable", sketchBlendEnable);
		sketch.setAttr("blend_mode",
				new Atom[] { Atom.newAtom(6), Atom.newAtom(7) });
		sketch.setAttr("antialias", sketchAntialias);
		sketch.setAttr(
				"glclearcolor",
				new Atom[] { Atom.newAtom(0.), Atom.newAtom(0.),
						Atom.newAtom(0.), Atom.newAtom(1.) });
		sketch.setAttr("fsaa", sketchFsaa);
		sketch.send("automatic", 0); /*
									 * set to not-automatic, to be able to use
									 * begin_capture and drawimmediate for
									 * capturing jit.gl.sketch as texture
									 */

		// instantiate Jitter texture object
		texture = new JitterObject("jit.gl.texture");
		texture.setAttr("drawto", context);
		texture.setAttr(
				"dim",
				new Atom[] { Atom.newAtom(texture_width),
						Atom.newAtom(texture_height) });
		// TODO: resize dimension

		// set arrays
		abPoint = new float[aNumMax][abPointCount][3];
		mbPoint = new float[aNumMax][abPointCount][3];
		abHeight = new float[aNumMax];
		abAngle = new float[aNumMax];
		rayCount = new int[aNumMax];
		mRayCount = new int[aNumMax];
		noiseCount = new int[aNumMax + 1]; // + 1 for the main noise field
		noiseF = new int[aNumMax + 1];
		noiseP = new int[aNumMax + 1];
		noiseStep = new int[aNumMax + 1];
		noiseWeight = new float[aNumMax + 1][3];
		rayHeight = new float[aNumMax];
		mRayHeight = new float[aNumMax];
		rayWidth = new float[aNumMax];
		rayColorOffset = new float[aNumMax];
		rayColorMult = new float[aNumMax];
		for (int c = 0; c < aNumMax; c++) {
			abHeight[c] = 0.5f;
			abAngle[c] = c * 0.2f;
			rayCount[c] = 50;
			rayHeight[c] = 1.f;
			rayWidth[c] = 1.f;
			rayColorOffset[c] = 0.4f;
			rayColorMult[c] = 0.5f;
		}
		for (int c = 0; c < aNumMax + 1; c++) {
			for (int xyz = 0; xyz < 2; xyz++) noiseWeight[c][xyz] = 0.f;
			noiseCount[c] = rayCount[0];
			noiseF[c] = 10;
			noiseP[c] = 0;
			noiseStep[c] = 1;
		}

		setMainBezier(); // sets outer bezier points to standard position
		setBezierPoints(); // define inner bezier points based on outer

		setNoise(); // set noise
	}

	/*
	 * draws and captures the aurora drawing to sketch object, and outputs
	 * sketch as jitter texture object
	 */
	public void bang() {
		texture.call("begin_capture"); // begin capturing	
		draw(); // draw aurora to sketch
		texture.call("end_capture"); // end capturing
		texture.call("draw"); // to output texture?
		outlet(0, "jit_gl_texture", texture.getAttr("name")); // output texture
	}

	/* draw aurora to jitter sketch object */
	public void draw() {

		if (noiseAnimation) moveNoise();

		
		sketch.call("reset");	// start drawing by resetting sketch object
		
		if (lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");

		for (int c = 0; c < aNum; c++) {
			// create local variables, to avoid conflict when life-updating
			// variables while rendering
			int _mbpC = abPointCount;
			float _mbP[][] = new float[_mbpC][3];
			for (int i = 0; i < _mbpC; i++) {
				for (int j = 0; j < 3; j++) {
					mbPoint[c][i][j] = (morphing) ? mbPoint[c][i][j] += (abPoint[c][i][j] - mbPoint[c][i][j])
							* slew
							: abPoint[c][i][j];
					_mbP[i][j] = mbPoint[c][i][j];
				}
			}
			mRayCount[c] = (morphing) ? (int) (mRayCount[c] += (rayCount[c] - mRayCount[c])
					* slew)
					: rayCount[c];
			int _r = mRayCount[c];
			mRayHeight[c] = (morphing) ? mRayHeight[c] += (rayHeight[c] - mRayHeight[c])
					* slew
					: rayHeight[c];
			float _rh = mRayHeight[c];
			float _rw = rayWidth[c];
			int _rs = raySegments;
			float _rco = rayColorOffset[c];
			float _rcm = rayColorMult[c];
			float _size[] = new float[3];
			float _pos[] = new float[3];
			float _vp[] = new float[3];
			for (int i = 0; i < 3; i++) {
				mSize[i] = (morphing) ? mSize[i] += (aSize[i] - mSize[i])
						* slew : aSize[i];
				_size[i] = mSize[i];
				mPos[i] = (morphing) ? mPos[i] += (aPos[i] - mPos[i]) * slew
						: aPos[i];
				_pos[i] = mPos[i];
				mVanishingPoint[i] = (morphing) ? mVanishingPoint[i] += (aVanishingPoint[i] - mVanishingPoint[i])
						* slew
						: aVanishingPoint[i];
				_vp[i] = mVanishingPoint[i];
			}
			
			Atom[] co = auroraColor(0.5f);
			sketch.call("glcolor", co);
			sketch.call("gllinewidth", _rw);

			float add = 1.f / (float) (_r - 1);

			for (int r = 0; r < _r; r++) {
				sketch.call("glbegin", "line_strip");

				float[] _p = P(add * r, _mbP); // returns array with xyz of main
												// point

				float _noisea = getNoise(c + 1, 0, r / (float) _r) * noiseWeight[c + 1][0];
				_noisea += getNoise(0, 0, r / (float) _r) * noiseWeight[0][0];		// add main noise field
				float _noiseb = getNoise(c + 1, 1, r / (float) _r) * noiseWeight[c + 1][1];
				_noiseb += getNoise(0, 1, r / (float) _r) * noiseWeight[0][1];		// add main noise field
				float _noisex = (float) (_noisea * Math.sin(aBeta) + _noiseb * Math.cos(aBeta));
				float _noisez = (float) (_noisea * Math.cos(aBeta) + _noiseb * Math.sin(aBeta));

				float _nx1 = _pos[0] + (_p[0] + _noisex) * _size[0];
				float _ny1 = _pos[1] + (_p[1]) * _size[1];
				float _nz1 = _pos[2] + (_p[2] + _noisez) * _size[2];

				float[] _n1 = { _nx1, _ny1, _nz1 };
				float[] _n2 = (_vp[1] > 0) ? vectorTowards(_n1, _vp, _rh
						* _size[1]) : vectorTowards(_n1, _vp, -_rh * _size[1]);

				for (int seg = 0; seg <= _rs; seg++) {
					float _segx = _nx1 + seg * (_n2[0] - _nx1)
							* (1.f / (float) _rs);
					float _segy = _ny1 + seg * (_n2[1] - _ny1)
							* (1.f / (float) _rs);
					float _segz = _nz1 + seg * (_n2[2] - _nz1)
							* (1.f / (float) _rs);

					sketch.call("glvertex", new Atom[] { Atom.newAtom(_segx),
							Atom.newAtom(_segy), Atom.newAtom(_segz) });
					// TODO: aSize seems to influence color
					Atom[] _segc = (rayColorByHeight) ? auroraColor((_segy
							- _rco - _pos[1])
							* _rcm) : auroraColor(seg / (float) _rs);
					sketch.call("glcolor", _segc);

				}

				sketch.call("glend");
			}
			
		}
		sketch.call("drawimmediate");	// call drawimmediate, to execute drawing of sketch object
	}

	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === === aurora functions = === === === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */

	/* set main bezier curve to standard symmetric position */
	void setMainBezier() {

		float l = 1.f / 2.0f;
		for (int c = 0; c < aNum; c++) {
			abPoint[c][0][0] = -l; // x 1st point
			abPoint[c][0][1] = 0.f; // y
			abPoint[c][0][2] = 0.f; // z

			abPoint[c][abPointCount - 1][0] = l; // x last point
			abPoint[c][abPointCount - 1][1] = 0.f; // y
			abPoint[c][abPointCount - 1][2] = 0.f; // z
		}

	}

	/* */
	void setBezierPoints() {

		for (int c = 0; c < aNumMax; c++) {

			// tip over abHeight by abAngle
			float newxz = (float) (Math.sin(abAngle[c]) * abHeight[c]);
			float newh = (float) (Math.cos(abAngle[c]) * abHeight[c]);

			// distance between first and last curve point
			float dx = abPoint[c][abPointCount - 2][0] - abPoint[c][0][0];
			float dy = abPoint[c][abPointCount - 2][1] - abPoint[c][0][1];
			float dz = abPoint[c][abPointCount - 2][2] - abPoint[c][0][2];
			float dl = (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
			float dbird = (float) Math.sqrt(dx * dx + dz * dz);

			aBeta = (float) Math.asin(dz / dbird);

			// define tipped over abHeight vector, by distance between endpoints
			float newx = newxz * dz / dl;
			float newz = newxz * dx / dl;

			// second bezier point
			abPoint[c][1][0] = abPoint[c][0][0] + newx;
			abPoint[c][1][1] = abPoint[c][0][1] + newh;
			abPoint[c][1][2] = abPoint[c][0][2] + newz;

			// if more than 2 inbetween bezier points

			// second to last bezier point
			abPoint[c][abPointCount - 2][0] = abPoint[c][abPointCount - 1][0] + newx;
			abPoint[c][abPointCount - 2][1] = abPoint[c][abPointCount - 1][1] + newh;
			abPoint[c][abPointCount - 2][2] = abPoint[c][abPointCount - 1][2] + newz;

		}

	}

	/* return color of aurora */
	private Atom[] auroraColor(float i) {
		float r = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getRed() / 255.f;
		float g = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getGreen() / 255.f;
		float b = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getBlue() / 255.f;
		return new Atom[] { Atom.newAtom(r), Atom.newAtom(g), Atom.newAtom(b) };
	}

	private float auroraSpectrum(float v) {
		float k = 0.5f - v * 0.5f; // 0.=red 0.5=green
		k = (k > 0.4f) ? k *= 0.75f : k; // TODO: finetune coloring, too choppy
											// right now
		k = (k < 0.25f && k > 0.125f) ? k * 1.5f : k;
		return restrict(k, 0.f, 0.5f);
		// float m = (k>1) ? 1 + ((float) Math.sin((k-1.f)*1.570796)) : 1 -
		// (float) Math.cos(k*1.570796);
		// return m/2.0f;
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

	private void setNoise() {
		
		int maxNoise = Math.max(noiseCount[0]+noiseF[0],Math.max(noiseCount[1]+noiseF[1],Math.max(noiseCount[2]+noiseF[2],noiseCount[3]+noiseF[3]))) + 1;
		noise = new float[aNumMax+1][2][maxNoise];

		// plus one, for main noise field, and all individual one
		for (int c = 0; c < aNumMax + 1; c++) {
			
//			noise[c][0] = new float[noiseCount[c] + noiseF[c] + 1];
//			noise[c][1] = new float[noiseCount[c] + noiseF[c] + 1];

			// random seed in noiseF interval, needs to include first and last
//			// element!
			for (int i = 0; i < noiseCount[c] + noiseF[c] + 1; i += noiseF[c]) {
				noise[c][0][i] = (float) Math.random() * 2 - 1.f;
				noise[c][1][i] = (float) Math.random() * 2 - 1.f;
			}

			// interpolation
			for (int xyz = 0; xyz < 2; xyz++) {
				for (int i = 0; i <= noiseCount[c]; i += noiseF[c]) {
					float a1 = noise[c][xyz][i];
					float a2 = noise[c][xyz][i + noiseF[c]];
					for (int j = i + 1; j < i + noiseF[c]; j++) {
						float x = (j - i) / (float) (noiseF[c] + 1);
						noise[c][xyz][j] = cosineInterpolate(a1, a2, x);
					}
				}
			}
		}

	}

	public void moveNoise() {

		for (int c = 0; c < aNum + 1; c++) {

			noiseP[c] += noiseStep[c]; // TODO: make noiseStep float
			if (noiseP[c] >= noiseF[c]) {

				while (noiseP[c] >= noiseF[c])
					noiseP[c] -= noiseF[c];
				// move hole array
				for (int xyz = 0; xyz < 2; xyz++) {
					for (int i = 0; i <= noiseCount[c]; i++) {
						noise[c][xyz][i] = noise[c][xyz][i + noiseF[c]];
					}
				}

				// new random seed at new end
				for (int xyz = 0; xyz < 2; xyz++) {
					noise[c][xyz][noiseCount[c] + noiseF[c]] = (float) Math
							.random() * 2 - 1.f;
				}

				// interpolate new segment
				for (int xyz = 0; xyz < 2; xyz++) {
					float a1 = noise[c][xyz][noiseCount[c]];
					float a2 = noise[c][xyz][noiseCount[c] + noiseF[c]];

					for (int j = noiseCount[c] + 1; j < noiseCount[c]
							+ noiseF[c]; j++) {
						float x = (j - noiseCount[c]) / (float) (noiseF[c] + 1);
						noise[c][xyz][j] = cosineInterpolate(a1, a2, x);
					}
				}
			}
		}
	}

	private float cosineInterpolate(float a, float b, float x) {
		float ft = x * 3.1415927f;
		float f = (1.f - (float) Math.cos(ft)) * .5f;
		return a * (1 - f) + b * f;
	}

	private float getNoise(int c, int xyz, float t) {
		int p = (int) (t * noiseCount[c]);
		// p = noiseCount - p;
		float m = 1.f;
		if (noiseSmoothEdge) {
			if (p < noiseEdge) {
				m = p / (float) noiseEdge;
			} else if (p > noiseCount[c] - noiseEdge) {
				m = (noiseCount[c] - p - 1) / (float) noiseEdge;
			}
		}
		try {
			return noise[c][xyz][p + noiseP[c]] * m;
		} catch(Exception e) {
			if(debug && p==0) post("error   c:"+c+"  xyz:"+xyz+"  p:"+p);
			return 0.f;
		}
		
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

	private float[] vectorTowards(float[] p1, float[] p2, float len) {
		float[] distance = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] }; // distance
																			// vector,
																			// p1
																			// to
																			// p2
		float distanceLength = (float) Math.sqrt(distance[0] * distance[0]
				+ distance[1] * distance[1] + distance[2] * distance[2]);
		float fact = len / distanceLength;
		distance[0] *= fact;
		distance[1] *= fact;
		distance[2] *= fact;
		float[] p = { p1[0] + distance[0], p1[1] + distance[1],
				p1[2] + distance[2] };
		return p;
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

	/* === === === === === === === aurora === === === === === === === === */

	/* set center position of rainbow */
	public void position(float x, float y, float z) {
		aPos[0] = x;
		aPos[1] = y;
		aPos[2] = z;
	}

	/* set main scale of rainbow */
	public void size(float x, float y, float z) {
		aSize[0] = x;
		aSize[1] = y;
		aSize[2] = z;
	}

	public void pA(float x, float y, float z) {
		for (int c = 0; c < aNum; c++) {
			abPoint[c][0][0] = x;
			abPoint[c][0][1] = y;
			abPoint[c][0][2] = z;
		}
		setBezierPoints();
	}

	public void pB(float x, float y, float z) {
		for (int c = 0; c < aNum; c++) {
			abPoint[c][abPointCount - 1][0] = x;
			abPoint[c][abPointCount - 1][1] = y;
			abPoint[c][abPointCount - 1][2] = z;
		}
		setBezierPoints();
	}

	public void raycount(int v) {
		raycount(v, v, v);
	}

	public void raycount(int v1, int v2, int v3) {
		rayCount[0] = (v1 > 2) ? v1 : 2;
		rayCount[1] = (v2 > 2) ? v2 : 2;
		rayCount[2] = (v3 > 2) ? v3 : 2;
		noiseCount[0] = rayCount[0];
		noiseCount[1] = rayCount[0];
		noiseCount[2] = rayCount[1];
		noiseCount[3] = rayCount[2];
		setNoise();
	}

	public void rayheight(float v) {
		rayheight(v, v, v);
	}

	public void rayheight(float v1, float v2, float v3) {
		rayHeight[0] = v1;
		rayHeight[1] = v2;
		rayHeight[2] = v3;
	}

	public void raywidth(float v) {
		raywidth(v, v, v);
	}

	public void raywidth(float v1, float v2, float v3) {
		rayWidth[0] = (v1 < 1) ? 1 : v1;
		rayWidth[1] = (v2 < 1) ? 1 : v2;
		rayWidth[2] = (v3 < 1) ? 1 : v3;
	}

	public void raysegments(int v) {
		raySegments = (v > 2) ? v : 2;
	}

	public void vanishingpoint(float x, float y, float z) {
		aVanishingPoint[0] = x;
		aVanishingPoint[1] = y;
		aVanishingPoint[2] = z;
	}

	// public void curveLength(float v) {
	// abLength = (v>0) ? v : 0.01f;
	// setMainBezier();
	// setBezierPoints();
	// }

	public void curveAngle(float v) {
		curveAngle(v, v, v);
	}

	public void curveAngle(float v1, float v2, float v3) {
		abAngle[0] = v1;
		abAngle[1] = v2;
		abAngle[2] = v3;
		setBezierPoints();
	}

	public void curveHeight(float v) {
		curveHeight(v, v, v);
	}

	public void curveHeight(float v1, float v2, float v3) {
		abHeight[0] = v1;
		abHeight[1] = v2;
		abHeight[2] = v3;
		setBezierPoints();
	}

	public void colorbyheight(int v) {
		rayColorByHeight = (v == 1) ? true : false;
	}

	public void coloroffset(float v) {
		coloroffset(v, v, v);
	}

	public void coloroffset(float v1, float v2, float v3) {
		rayColorOffset[0] = v1;
		rayColorOffset[1] = v2;
		rayColorOffset[2] = v3;
	}

	public void colormult(float v) {
		colormult(v, v, v);
	}

	public void colormult(float v1, float v2, float v3) {
		rayColorMult[0] = v1;
		rayColorMult[1] = v2;
		rayColorMult[2] = v3;
	}

	/* toggle morphing */
	public void morph(int v) {
		morphing = (v == 1) ? true : false;
	}

	/* set morphing speed */
	public void morphspeed(float v) {
		slew = (v > 0) ? v : 0.01f;
	}

	public void noisef(int v) {
		noisef(v, v, v, v);
	}

	public void noisef(int v0, int v1, int v2, int v3) {
		noiseF[0] = (v0 > 1) ? v0 : 1;
		noiseF[1] = (v1 > 1) ? v1 : 1;
		noiseF[2] = (v2 > 1) ? v2 : 1;
		noiseF[3] = (v3 > 1) ? v3 : 1;
		setNoise();
	}
	
	public void noiseweight(float v0, float v1, float v2, float v3) {
		noiseweight(v0,v0,v1,v1,v2,v2,v3,v3);
	}

	public void noiseweight(float vx0, float vz0, float vx1, float vz1, float vx2,float vz2, float vx3, float vz3) {
		noiseWeight[0][0] = vx0;
		noiseWeight[0][1] = vz0;
		noiseWeight[1][0] = vx1;
		noiseWeight[1][1] = vz1;
		noiseWeight[2][0] = vx2;
		noiseWeight[2][1] = vz2;
		noiseWeight[3][0] = vx3;
		noiseWeight[3][1] = vz3;
	}

	public void noiseanimation(int v) {
		noiseAnimation = (v == 1) ? true : false;
	}

	public void noiseedge(int v) {
		noiseSmoothEdge = (v <= 0) ? false : true;
		noiseEdge = (v > 0) ? v : 0;
	}

	public void noisestep(int v) {
		noisestep(v, v, v, v);
	}

	public void noisestep(int v0, int v1, int v2, int v3) {
		noiseStep[0] = (v0 > 0) ? v0 : 1;
		noiseStep[1] = (v1 > 0) ? v1 : 1;
		noiseStep[2] = (v2 > 0) ? v2 : 1;
		noiseStep[3] = (v3 > 0) ? v3 : 1;
	}
	
	public void num(int v) {
		aNum = (v<aNumMax) ? (v>0) ? v : 1 : aNumMax;
	}

}
