/*
 * Figure de la terre - mxj objects
 * 
 * aurora borealis
 * main bezier curve
 * 
 */

import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;


public class FigureAurora extends MaxObject {
	
	
	// render context for jitter opengl objects
	private String context = "foo";
	
	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	private int texture_width = 640;	// TODO: include resize function
	private int texture_height = 480;
	
	// Bezier curve, aurora backbone
	private float[] aPos = {0,0,0};
	private float aBeta = 0.f;			// birdview angle between outer bezier points
	private float[] aSize = {1,1,1};
	private int abPointCount = 4;		// main bezier point count, leave count flexible just in case
	private float abPoint[][];			// main bezier points
	private float abLength = 0.5f;		// distance between outer points
	private float abHeight = 0.5f;		// z of inner bezier points defines height of curve
	private float abAngle = 0.0f;		// angle of curve, leaning right or left ?
	
	// aurora rays
	private int rayCount = 100;
	private float rayHeight = 1.f;
	private float rayWidth = 1.f;
	private int raySegments = 4;
	private boolean rayColorByHeight = true;
	private float rayColorOffset = 0.4f;
	private float rayColorMult = 0.5f;
	private float[] aVanishingPoint = {0,100,0};
	
	// noise
	private int noiseCount = 100;
	private int noiseF = 10;	// noise frequency
	private int noiseP = 0;		// noise pointer, for animating the noise
	private int noiseStep = 2;
	private boolean noiseAnimation = true;
	private int noiseEdge = 0;
	private boolean noiseSmoothEdge = false;
	private float[] noisea;
	private float[] noiseb;
	private float noiseWeight[] = {0.1f,0.01f};
	
	// morphing
	private boolean morphing = false;
	private float slew = 0.1f;			// morphing speed
	private int mRayCount = 100;
	private float mPos[] = {0,0,0};
	private float mSize[] = {1,1,1};
	private float mbPoint[][];
	private float[] mVanishingPoint = {0,100,0};
//	private float mbHeight = 0.5f;
//	private float mbAngle = 0.0f;
	private float mRayHeight = 1.f;
	
	// display parameters
	private boolean lineSmooth = false;
	private int sketchDepthEnable = 0;
	private int sketchAntialias = 0;
	private int sketchFsaa = 0;
	private int sketchBlendEnable = 0;

	
	/* === === === === === === === main functions === === === === === === === === */
	
	
	
	/* instantiating mxj without argument */
	public FigureAurora() {
		bail("please provide render context as argument");
	}
	
	/* instantiate mxj with render context as argument */
	public FigureAurora(String c) {
		
		context = c;		// render context
		declareIO(2,1);		// declare 2 inlets, 1 outlet of DataTypes.ALL
		
		// assist message for inlets and outlets (mouse hover)
		setInletAssist(new String[] {"bang to compute and draw", "input settings"});
		setOutletAssist(new String[] {"outputs jit.gl.texture object"});
		
		// instantiate Jitter sketch object
		sketch = new JitterObject("jit.gl.sketch");
		sketch.setAttr("drawto", context);
		sketch.setAttr("depth_enable", sketchDepthEnable);
		sketch.setAttr("blend_enable", sketchBlendEnable);
		sketch.setAttr("blend_mode", new Atom[]{Atom.newAtom(6), Atom.newAtom(7)});
		sketch.setAttr("antialias",sketchAntialias);
		sketch.setAttr("glclearcolor", new Atom[]{Atom.newAtom(0.),Atom.newAtom(0.),Atom.newAtom(0.),Atom.newAtom(1.)});
		sketch.setAttr("fsaa", sketchFsaa);
		sketch.send("automatic", 0);	/* set to not-automatic, to be able to use
										   begin_capture and drawimmediate 
										   for capturing jit.gl.sketch as texture */
		
		
		// instantiate Jitter texture object
		texture = new JitterObject("jit.gl.texture");
		texture.setAttr("drawto", context);
		texture.setAttr("dim", new Atom[]{Atom.newAtom(texture_width),Atom.newAtom(texture_height)}); 
		// TODO: resize dimension

		// set arrays
		abPoint = new float[abPointCount][3];
		mbPoint = new float[abPointCount][3];
		
		setMainBezier();	// sets outer bezier points to standard position
		setBezierPoints();	// define inner bezier points based on outer
		
		setNoise();			// set noise
	}
	
	/* draws and captures the aurora drawing to sketch object, 
	 * and outputs sketch as jitter texture object */
	public void bang() {
		texture.call("begin_capture");			// begin capturing	
		draw();									// draw aurora to sketch
		texture.call("end_capture");			// end capturing
		texture.call("draw");					// to output texture? 
		outlet(0,"jit_gl_texture",texture.getAttr("name"));		// output texture
	}
	
	/* draw aurora to jitter sketch object */
	public void draw() {
		
		if(noiseAnimation) moveNoise();
		
		// create local variables, to avoid conflict when life-updating variables while rendering
		int _mbpC = abPointCount;
		float _mbP[][] = new float[_mbpC][3];
		for(int i=0; i<_mbpC; i++) {
			for(int j=0; j<3; j++) {
				mbPoint[i][j] = (morphing) ? mbPoint[i][j] += (abPoint[i][j]-mbPoint[i][j])*slew : abPoint[i][j];
				_mbP[i][j] = mbPoint[i][j];
			}
		}
		mRayCount = (morphing) ? (int) (mRayCount += (rayCount-mRayCount)*slew) : rayCount;
		int _r = mRayCount;
		mRayHeight = (morphing) ? mRayHeight += (rayHeight-mRayHeight)*slew : rayHeight;
		float _rh = mRayHeight;
		float _rw = rayWidth;
		int _rs = raySegments;
		float _rco = rayColorOffset;
		float _rcm = rayColorMult;
		
		float _size[] = new float[3];
		float _pos[] = new float[3];
		float _vp[] = new float[3];
		for(int i=0; i<3; i++) {
			mSize[i] = (morphing) ? mSize[i] += (aSize[i]-mSize[i])*slew : aSize[i];
			_size[i] = mSize[i];
			mPos[i] = (morphing) ? mPos[i] += (aPos[i]-mPos[i])*slew : aPos[i];
			_pos[i] = mPos[i];
			mVanishingPoint[i] = (morphing) ? mVanishingPoint[i] += (aVanishingPoint[i]-mVanishingPoint[i])*slew : aVanishingPoint[i];
			_vp[i] = mVanishingPoint[i];
		}
		
		// start drawing by reseting sketch object
		sketch.call("reset");
		
		if(lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");
		
		Atom[] c = auroraColor(0.5f);
		sketch.call("glcolor", c);
		sketch.call("gllinewidth", _rw);
		
		
		float add = 1.f / (float) (_r-1);
		
		for(int r=0; r<_r; r++) {
			sketch.call("glbegin", "line_strip");
			
			float[] _p = P(add*r, _mbP);			// returns array with xyz of main point
			
			float _noisea = getNoiseA(r/(float)_r)*noiseWeight[0];
			float _noiseb = getNoiseB(r/(float)_r)*noiseWeight[1];
			float _noisex = (float) (_noisea*Math.sin(aBeta) + _noiseb*Math.cos(aBeta));
			float _noisez = (float) (_noisea*Math.cos(aBeta) + _noiseb*Math.sin(aBeta));
			
			float _nx1 = _pos[0] + (_p[0]+_noisex)*_size[0];
			float _ny1 = _pos[1] + (_p[1])*_size[1];
			float _nz1 = _pos[2] + (_p[2]+_noisez)*_size[2];
			
			float[] _n1 = { _nx1, _ny1, _nz1 };
			float[] _n2 = (_vp[1]>0) ? vectorTowards(_n1, _vp, _rh*_size[1]) : vectorTowards(_n1, _vp, -_rh*_size[1]);
			
			for(int seg=0; seg<=_rs; seg++) {
				float _segx = _nx1 + seg*(_n2[0]-_nx1)*(1.f/(float) _rs);
				float _segy = _ny1 + seg*(_n2[1]-_ny1)*(1.f/(float) _rs);
				float _segz = _nz1 + seg*(_n2[2]-_nz1)*(1.f/(float) _rs);
				
				sketch.call("glvertex", new Atom[]{Atom.newAtom(_segx),Atom.newAtom(_segy),Atom.newAtom(_segz)});
				// TODO: aSize seems to influence color
				Atom[] _segc = (rayColorByHeight) ? auroraColor((_segy-_rco-_pos[1])*_rcm) : auroraColor(seg/(float) _rs);
				sketch.call("glcolor", _segc);
				
				
			}
			
			
			sketch.call("glend");
		}
		
		
		// call drawimmediate, to execute drawing of sketch object
		sketch.call("drawimmediate");		
		
		
	}
	
	
	
	
	/* === === === === === === === === ====== === === === === === === === === === */
	/* === === === === === === =  aurora functions  = === === === === === === === */
	/* === === === === === === === === ====== === === === === === === === === === */
	
	
	/* set main bezier curve to standard symmetric position */
	void setMainBezier() {
		
		float l = abLength/2.0f;
		
		abPoint[0][0] = -l;	// x 	1st point
		abPoint[0][1] = 0.f;	// y
		abPoint[0][2] = 0.f;	// z
		
		abPoint[abPointCount-1][0] = l;		// x	last point
		abPoint[abPointCount-1][1] = 0.f;	// y
		abPoint[abPointCount-1][2] = 0.f;	// z

	}
	
	/* */
	void setBezierPoints() {
		
		
		// tip over abHeight by abAngle
		float newxz = (float) (Math.sin(abAngle)*abHeight);
		float newh = (float) (Math.cos(abAngle)*abHeight);
		
		// distance between first and last curve point
		float dx = abPoint[abPointCount-2][0] - abPoint[0][0];
		float dy = abPoint[abPointCount-2][1] - abPoint[0][1];
		float dz = abPoint[abPointCount-2][2] - abPoint[0][2];
		float dl = (float) Math.sqrt(dx*dx + dy*dy + dz*dz);
		float dbird = (float) Math.sqrt(dx*dx + dz*dz);
		
		aBeta = (float) Math.asin(dz/dbird);
		
		// define tipped over abHeight vector, by distance between endpoints
		float newx = newxz * dz/dl;
		float newz = newxz * dx/dl;
		
		
		
		
		// second bezier point
		abPoint[1][0] = abPoint[0][0] + newx;
		abPoint[1][1] = abPoint[0][1] + newh;
		abPoint[1][2] = abPoint[0][2] + newz;
		
		// if more than 2 inbetween bezier points
		
		// second to last bezier point
		abPoint[abPointCount-2][0] = abPoint[abPointCount-1][0] + newx;
		abPoint[abPointCount-2][1] = abPoint[abPointCount-1][1] + newh;
		abPoint[abPointCount-2][2] = abPoint[abPointCount-1][2] + newz;
		
	}
	
	
	/* return color of aurora  */
	private Atom[] auroraColor(float i) {
		float r = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
		float g = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
		float b = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
		return new Atom[]{Atom.newAtom(r),Atom.newAtom(g),Atom.newAtom(b)};
	}
	
	private float auroraSpectrum(float v) {
		float k = 0.5f - v*0.5f;					// 0.=red    0.5=green
		k = (k>0.4f) ? k*=0.75f : k;				// TODO: finetune coloring, too choppy right now
		k = (k<0.25f && k > 0.125f) ? k*1.5f : k;
		return restrict(k,0.f,0.5f);
//		float m = (k>1) ? 1 + ((float) Math.sin((k-1.f)*1.570796)) : 1 - (float) Math.cos(k*1.570796);
//		return m/2.0f;
	}
	
	
	/* === === === === === === === === ====== === === === === === === === === === */
	/* === === === === === ===    bezier functions    === === === === === === === */
	/* === === === === === === === === ====== === === === === === === === === === */
	
	
	
	
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
	
	
	
	private void setNoise() {
		
		noiseP = 0;
		noisea = new float[noiseCount+noiseF+1];	// noiseC is multiple of noiseF   o--o--o--o
		noiseb = new float[noiseCount+noiseF+1];
		
		// random seed in noiesF interval, needs to include first and last element!
		for(int i=0; i<noiseCount+noiseF+1; i+=noiseF) {
			noisea[i] = (float) Math.random()*2 - 1.f;
			noiseb[i] = (float) Math.random()*2 - 1.f;
		}
		
		// interpolation
		for(int i=0; i<=noiseCount; i+=noiseF) {
			float a1 = noisea[i];
			float a2 = noisea[i+noiseF];
			float b1 = noiseb[i];
			float b2 = noiseb[i+noiseF];
			for(int j=i+1; j<i+noiseF; j++) {
				float x = (j-i)/(float) (noiseF+1);
				noisea[j] = cosineInterpolate(a1,a2,x);
				noiseb[j] = cosineInterpolate(b1,b2,x);
			}
		}
		
		// additional smoothing  TODO
		
		
		// smooth edges   TODO
		
		
	}
	
	public void moveNoise() {
		noiseP+=noiseStep;
		if(noiseP>=noiseF) {
			
			while(noiseP>=noiseF) noiseP-=noiseF;
			//  move hole array
			for(int i=0; i<=noiseCount; i++) {
				noisea[i] = noisea[i+noiseF];
				noiseb[i] = noiseb[i+noiseF];
			}
			
			// new random seed at new end
			noisea[noiseCount+noiseF] = (float) Math.random()*2 - 1.f;
			noiseb[noiseCount+noiseF] = (float) Math.random()*2 - 1.f;
			
			// interpolate new segment
			float a1 = noisea[noiseCount];
			float a2 = noisea[noiseCount+noiseF];
			float b1 = noiseb[noiseCount];
			float b2 = noiseb[noiseCount+noiseF];
			for(int j=noiseCount+1; j<noiseCount+noiseF; j++) {
				float x = (j-noiseCount)/(float) (noiseF+1);
				noisea[j] = cosineInterpolate(a1,a2,x);
				noiseb[j] = cosineInterpolate(b1,b2,x);
			}

		}
	}
	
	private float cosineInterpolate(float a, float b, float x) {
		float ft = x*3.1415927f;
		float f = (1.f - (float) Math.cos(ft)) * .5f;
		return a*(1-f) + b*f;
	}
	
	private float getNoiseA(float t) {
		int p = (int) (t*noiseCount);
//		p = noiseCount - p;
		float m = 1.f;
		if(noiseSmoothEdge) {
			if(p<noiseEdge) {
				m = p / (float) noiseEdge;
			} else if(p > noiseCount-noiseEdge) {
				m = (noiseCount-p-1) / (float) noiseEdge;
			}
		}
		return noisea[p + noiseP] * m;
	}
	
	private float getNoiseB(float t) {
		int p = (int) (t*noiseCount);
//		p = noiseCount - p;
		float m = 1.f;
		if(noiseSmoothEdge) {
			if(p<noiseEdge) {
				m = p / (float) noiseEdge;
			} else if(p > noiseCount-noiseEdge) {
				m = (noiseCount-p-1) / (float) noiseEdge;
			}
		}
		return noiseb[p + noiseP] * m;
	}
	
	
	
	/* === === === === === === === === ====== === === === === === === === === === */
	/* === === === === === == helpful little functions == === === === === === === */
	/* === === === === === === === === ====== === === === === === === === === === */
	
	private float restrict(float v) {
		return restrict(v,-1.f,1.f);
	}
	
	private float restrict(float v, float min, float max) {
		if(v<min) return min;
		if(v>max) return max;
		return v;
	}
	
	private float[] vectorTowards( float[] p1, float[] p2, float len) {
		float[] distance = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };	// distance vector, p1 to p2
		float distanceLength = (float) Math.sqrt(distance[0]*distance[0] + distance[1]*distance[1] + distance[2]*distance[2]);
		float fact = len/distanceLength;
		distance[0]*=fact;
		distance[1]*=fact;
		distance[2]*=fact;
		float[] p = { p1[0]+distance[0], p1[1]+distance[1], p1[2]+distance[2] };
		return p;
	}
	
	
	
	/* === === === === === === === === ====== === === === === === === === === === */
	/* === === === === === === === set parameters === === === === === === === === */
	/* === === === === === === === === ====== === === === === === === === === === */
	
	
	
	
	
	/* === === === === === === === jit.gl.sketch  === === === === === === === === */
	
	public void depthEnable(int v) {
		sketchDepthEnable = (v==1) ? 1 : 0;
		sketch.setAttr("depth_enable", sketchDepthEnable);
	}
	
	public void blendEnable(int v) {
		sketchBlendEnable = (v==1) ? 1 : 0;
		sketch.setAttr("blend_enable", sketchBlendEnable);
	}
	
	public void blendMode(int b1, int b2) {
		sketch.setAttr("blend_mode", new Atom[]{Atom.newAtom(b1), Atom.newAtom(b2)});
	}
	
	public void antialias(int v) {
		sketchAntialias = (v==1) ? 1 : 0;
		sketch.setAttr("antialias", sketchAntialias);
	}
	
	public void clearColor(float r, float g, float b, float a) {
		sketch.setAttr("glclearcolor", new Atom[]{Atom.newAtom(r),Atom.newAtom(g),Atom.newAtom(b),Atom.newAtom(a)});
	}
	
	public void fsaa(int v) {
		sketchFsaa = (v==1) ? 1 : 0;
		sketch.setAttr("fsaa", sketchFsaa);
	}
	
	public void linesmooth(int v) {
	    lineSmooth = (v==1) ? true : false;
	}
	
	
	
	/* === === === === === === ===    aurora    === === === === === === === === */
	
	
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
		abPoint[0][0] = x;
		abPoint[0][1] = y;
		abPoint[0][2] = z;
		setBezierPoints();
	}
	
	public void pB(float x, float y, float z) {
		abPoint[abPointCount-1][0] = x;
		abPoint[abPointCount-1][1] = y;
		abPoint[abPointCount-1][2] = z;
		setBezierPoints();
	}
	
	
	public void raycount(int v) {
		rayCount = (v>2) ? v : 2;
		noiseCount = rayCount;
		setNoise();
	}
	
	public void rayheight(float v) {
		rayHeight = v;
	}
	
	public void raywidth(float v) {
		rayWidth = (v<1) ? 1 : v;
	}
	
	public void raysegments(int v) {
		raySegments = (v>2) ? v : 2;
	}
	
	
	public void vanishingpoint(float x, float y, float z) {
		aVanishingPoint[0] = x;
		aVanishingPoint[1] = y;
		aVanishingPoint[2] = z;
	}
	
	
	public void curveLength(float v) {
		abLength = (v>0) ? v : 0.01f;
		setMainBezier();	
		setBezierPoints();
	}
	
	public void curveAngle(float v) {
		abAngle = v;
		setBezierPoints();
	}
	
	public void curveHeight(float v) {
		abHeight = v;
		setBezierPoints();
	}
	
	
	public void colorbyheight(int v) {
		rayColorByHeight = (v==1) ? true : false;
	}
	
	public void coloroffset(float v) {
		rayColorOffset = v;
	}
	
	public void colormult(float v) {
		rayColorMult = v;
	}
	
	/* toggle morphing */
	public void morph(int v) {
		morphing = (v==1) ? true : false;
	}
	
	/* set morphing speed */
	public void morphspeed(float v) {
		slew = (v>0) ? v : 0.01f;
	}
	
	
	public void noisef(int v) {
		noiseF = (v>1) ? v : 1;
		setNoise();
	}
	
	public void noiseweight(float v1, float v2) {
		noiseWeight[0] = v1;
		noiseWeight[1] = v2;
	}
	
	public void noiseanimation(int v) {
		noiseAnimation = (v==1) ? true : false;
	}
	
	public void noiseedge(int v) {
		noiseSmoothEdge = (v<=0) ? false : true;
		noiseEdge = (v>0) ? v : 0;
	}
	
	public void noisestep(int v) {
		noiseStep = (v>0) ? v : 1;
	}

}
