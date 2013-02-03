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
	private float[] aSize = {1,1,1};
	private int mbPointCount = 4;		// main bezier point count, leave count flexible just in case
	private float mbPoint[][];			// main bezier points
	private float mbLength = 0.5f;		// distance between outer points
	private float mbHeight = 0.5f;		// z of inner bezier points defines height of curve
	private float mbAngle = 0.1f;		// angle of curve, leaning right or left ?
	
	// aurora rays
	private int rayCount = 50;
	private float rayHeight = 1.f;
	private float rayWidth = 1.f;
	
	// display parameters
	private boolean lineSmooth = false;
	private int sketchDepthEnable = 0;
	private int sketchAntialias = 0;
	private int sketchFsaa = 0;

	
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
		mbPoint = new float[mbPointCount][3];
		
		setMainBezier();	// sets outer bezier points to standard position
		setBezierPoints();	// define inner bezier points based on outer
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
		
		// create local variables, to avoid conflict when life-updating variables while rendering
		int _mbpC = mbPointCount;
		float _mbP[][] = new float[_mbpC][3];
		for(int i=0; i<_mbpC; i++) {
			for(int j=0; j<3; j++) _mbP[i][j] = mbPoint[i][j];
		}
		int _r = rayCount;
		float _rh = rayHeight;
		float _rw = rayWidth;
		float _size[] = { aSize[0], aSize[1], aSize[2] };
		float _pos[] = { aPos[0], aPos[1], aPos[2] };
		
		
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
			
			float _nx = _pos[0] + _p[0]*_size[0];
			float _ny = _pos[1] + _p[1]*_size[1];
			float _nz = _pos[2] + _p[2]*_size[2];
			sketch.call("glvertex", new Atom[]{Atom.newAtom(_nx),Atom.newAtom(_ny),Atom.newAtom(_nz)});
			_ny = _pos[1] + (_p[1]+_rh)*_size[1];
			sketch.call("glvertex", new Atom[]{Atom.newAtom(_nx),Atom.newAtom(_ny),Atom.newAtom(_nz)});
			
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
		
		float l = mbLength/2.0f;
		
		mbPoint[0][0] = -l;	// x 	1st point
		mbPoint[0][1] = 0.f;	// y
		mbPoint[0][2] = 0.f;	// z
		
		mbPoint[mbPointCount-1][0] = l;		// x	last point
		mbPoint[mbPointCount-1][1] = 0.f;	// y
		mbPoint[mbPointCount-1][2] = 0.f;	// z

	}
	
	/* */
	void setBezierPoints() {
		
		
		// tip over mbHeight by mbAngle
		float newxz = (float) (Math.sin(mbAngle)*mbHeight);
		float newh = (float) (Math.cos(mbAngle)*mbHeight);
		
		// distance between first and last curve point
		float dx = mbPoint[mbPointCount-2][0] - mbPoint[0][0];
		float dy = mbPoint[mbPointCount-2][1] - mbPoint[0][1];
		float dz = mbPoint[mbPointCount-2][2] - mbPoint[0][2];
		float dl = (float) Math.sqrt(dx*dx + dy*dy + dz*dz);
		
		// define tipped over mbHeight vector, by distance between endpoints
		float newx = newxz * dz/dl;
		float newz = newxz * dx/dl;
		
		
		
		
		// second bezier point
		mbPoint[1][0] = mbPoint[0][0] + newx;
		mbPoint[1][1] = mbPoint[0][1] + newh;
		mbPoint[1][2] = mbPoint[0][2] + newz;
		
		// if more than 2 inbetween bezier points
		
		// second to last bezier point
		mbPoint[mbPointCount-2][0] = mbPoint[mbPointCount-1][0] + newx;
		mbPoint[mbPointCount-2][1] = mbPoint[mbPointCount-1][1] + newh;
		mbPoint[mbPointCount-2][2] = mbPoint[mbPointCount-1][2] + newz;
		
	}
	
	
	/* return color of aurora  */
	private Atom[] auroraColor(float i) {
		float r = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getRed()/255.f;
		float g = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getGreen()/255.f;
		float b = Color.getHSBColor(auroraSpectrum(i), 1.0f, 0.8f).getBlue()/255.f;
		return new Atom[]{Atom.newAtom(r),Atom.newAtom(g),Atom.newAtom(b)};
	}
	
	private float auroraSpectrum(float v) {
		float k = v*2.f;
		float m = (k>1) ? 1 + ((float) Math.sin((k-1.f)*1.570796)) : 1 - (float) Math.cos(k*1.570796);
		return m/2.0f;
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
	
	
	
	/* === === === === === === === === ====== === === === === === === === === === */
	/* === === === === === === === set parameters === === === === === === === === */
	/* === === === === === === === === ====== === === === === === === === === === */
	
	
	
	
	
	/* === === === === === === === jit.gl.sketch  === === === === === === === === */
	
	public void depthEnable(int v) {
		sketchDepthEnable = (v==1) ? 1 : 0;
		sketch.setAttr("depth_enable", sketchDepthEnable);
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
		mbPoint[0][0] = x;
		mbPoint[0][1] = y;
		mbPoint[0][2] = z;
		setBezierPoints();
	}
	
	public void pB(float x, float y, float z) {
		mbPoint[mbPointCount-1][0] = x;
		mbPoint[mbPointCount-1][1] = y;
		mbPoint[mbPointCount-1][2] = z;
		setBezierPoints();
	}
	
	
	public void raycount(int v) {
		rayCount = (v>2) ? v : 2;
	}
	
	public void rayheight(float v) {
		rayHeight = v;
	}
	
	public void raywidth(float v) {
		rayWidth = (v<1) ? 1 : v;
	}
	
	public void curveLength(float v) {
		mbLength = (v>0) ? v : 0.01f;
		setMainBezier();	
		setBezierPoints();
	}
	
	public void curveAngle(float v) {
		mbAngle = v;
		setBezierPoints();
	}
	
	public void curveHeight(float v) {
		mbHeight = v;
		setBezierPoints();
	}
	
	
	
	
	

}
