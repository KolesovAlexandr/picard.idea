package picard.util;

/**
 * Created by Aleksandr_Kolesov on 10/8/2015.
 */
public class LoopArray {

    private int _length;
    private int[] _array;
    private int[] _pass;
    public LoopArray(int length){
        _length = length;
        _array = new int[_length];
        _pass = new int[_length];

    }
    public int getIndex(int i){
        int index = i%_length;
        int pass = i/_length;
        if (_pass[index] != pass){
            _array[index]=0;
            _pass[index] = pass;
        }

        return index;
    }
    public void incriment(int i){
        _array[i]++;
    }
    public  int get(int i){
        return _array[i];
    }
    public  int lengt(){
        return _length;
    }
}
