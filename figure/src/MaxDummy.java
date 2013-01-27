import com.cycling74.max.*;

public class MaxDummy extends MaxObject {
	
	private int fun=0;
	
	public MaxDummy() {
//		declareInlets(new int[]{DataTypes.ANYTHING, DataTypes.ANYTHING});
		declareIO(2,1);
		declareAttribute("fun");
		post("MaxDummy class, zum Dienst!");
		error("Nicht alles ist gut.");
	}
	
	public void bang() {
		if(fun==0) {
			outlet(0, "Twat!");
		} else {
			outlet(0, "Twitteltuueee!");
		}
	}	
	
	public void inlet(int i) {
		if(getInlet() == 0) {
			outlet(0, "Eins");
		} else {
			outlet(0, "Zwei");
		}
	}
	
	private void mystery() {
		outlet(0, "Wuck!");
	}
}
