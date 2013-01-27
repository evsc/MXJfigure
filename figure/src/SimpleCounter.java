import com.cycling74.max.*;

public class SimpleCounter extends MaxObject implements Executable {

	private MaxClock clock;
	private int count=0;
	private float interval=100.f;
	
	public SimpleCounter() {
		declareAttribute("count");
		declareAttribute("interval");
		clock = new MaxClock(this);
	}
	
	public void execute() {
		bang();
		clock.delay(interval);
	}
	
	public void bang() {
		count++;
		outlet(0, count);
	}
	
	protected void notifyDeleted() {
		clock.unset();
	}
	
	public void inlet(int i) {
		if(i==1)
			clock.delay(0);
		else
			clock.unset();
	}
	
}
