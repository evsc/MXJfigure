/*
 * Figure de la terre - mxj objects
 * 
 * paisley wallpaper
 * simple grid
 * 
 */

import java.awt.Color;

import com.cycling74.max.*;
import com.cycling74.jitter.*;



public class FigurePaisley extends MaxObject {
	
	// render context for jitter opengl objects
	private String context = "foo";
	private boolean debug = true;

	// jitter objects
	JitterObject sketch;
	JitterObject texture;
	private int texture_width = 640; 		// TODO: include resize function
	private int texture_height = 480;
	
	
	
	// grid
	private float gridx = 0.2f;
	private float gridy = 0.25f;
	
	private float[] pPos = {0,0,0};	// center position
	private float[] pSize = {1,1,1};		// scale
	
	// 
	private int pNumMax = 100;	// maximum number of paisley patterns
	
	
	
	// morphing
	private boolean morphing = false;
	private float slew = 0.1f; 				// morphing speed
	private float[] mPos = { 0, 0, 0 };
	private float[] mSize = { 1, 1, 1 };
	
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
				new Atom[] { Atom.newAtom(0.), Atom.newAtom(0.), Atom.newAtom(0.), Atom.newAtom(1.) });
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

		sketch.call("reset");
		
		if (lineSmooth) sketch.call("glenable", "line_smooth");
		else sketch.call("gldisable", "line_smooth");
		
		// create local variables, to avoid conflict when life-updating
		// variables while rendering
		
		float _size[] = new float[2];		// scale from origin
		float _pos[] = new float[2];		// origin position
		
		for (int i = 0; i < 2; i++) {		// only need x and y values
			mSize[i] = (morphing) ? mSize[i] += (pSize[i] - mSize[i]) * slew : pSize[i];
			_size[i] = mSize[i];
			mPos[i] = (morphing) ? mPos[i] += (pPos[i] - mPos[i]) * slew : pPos[i];
			_pos[i] = mPos[i];
		}
		
		float _gx = gridx;
		float _gy = gridy;
		
		
		
		
		
		
		// - - - - - - - - - - - - - - - - - - - - - - -
		
		int _i = 0;			// current item
		boolean goon = true;	// go on
		
		float _x = (float) Math.floor(1.f / _gx) * -_gx;
		float _y = (float) Math.floor(0.7f / _gy) * -_gy;
		float borderx = -_x;
		float bordery = -_y;
		int row = 0;
		
		while(goon) {

			drawPaisley(_pos[0] + _x*_size[0],_pos[1] + _y*_size[1]);
			
			
			// adapt the x and y values to the paisley grid
			// every second row the elements are positioned in between the top row grid
			_x+= _gx*2;
			if(_x > borderx) {
				row++;
				_x = (row%2 == 0) ? -borderx : -borderx + _gx;
				_y += _gy;
				if(_y > bordery) goon = false;
			}

			_i++;
			if(_i > pNumMax) goon = false;
		}

		
		// call drawimmediate, to execute drawing of sketch object
		sketch.call("drawimmediate");	
		
	}
	
	
	public void drawPaisley(float _x, float _y) {
//		if(debug) post("drawPaisly( "+_x+", "+_y+")");
		
		float w = 0.1f;
		
		Atom[] bandColor = new Atom[]{Atom.newAtom(255),Atom.newAtom(1),Atom.newAtom(1)};
		sketch.call("glcolor", bandColor);
		
		sketch.call("glbegin", "polygon");
		
		sketch.call("glvertex", new Atom[]{	Atom.newAtom(_x+w), Atom.newAtom(_y+w), Atom.newAtom(0)});
		sketch.call("glvertex", new Atom[]{	Atom.newAtom(_x-w), Atom.newAtom(_y+w), Atom.newAtom(0)});
		sketch.call("glvertex", new Atom[]{	Atom.newAtom(_x-w), Atom.newAtom(_y-w), Atom.newAtom(0)});
		sketch.call("glvertex", new Atom[]{	Atom.newAtom(_x+w), Atom.newAtom(_y-w), Atom.newAtom(0)});
		
		sketch.call("glend");
		
	}

	
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === === paisley functions ==== === === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	
	
	
	
	
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
	
	
	
	
	
	
	/*
	 * === === === === === === === === ====== === === === === === === === ===
	 * === === === === === == helpful little functions == === === === === ===
	 * === === === === === === === === ====== === === === === === === === ===
	 */
	
	
	
	
	
	
	
	
	
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
		gridx = x;
		gridy = y;
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
		
		
		outlet(1,"script send gui_sizex set "+pSize[0]);
		outlet(1,"script send gui_sizey set "+pSize[1]);
		outlet(1,"script send gui_sizez set "+pSize[2]);
		outlet(1,"script send gui_positionx set "+pPos[0]);
		outlet(1,"script send gui_positiony set "+pPos[1]);
		outlet(1,"script send gui_positionz set "+pPos[2]);
		
		outlet(1,"script send gui_gridx set "+gridx);
		outlet(1,"script send gui_gridy set "+gridy);
		outlet(1, "script send gui_grid set grid " + gridx + " " + gridy);
		
	}
	

	
	
	
	
}
