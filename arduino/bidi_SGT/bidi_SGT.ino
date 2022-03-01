void setup() {
  // put your setup code here, to run once:
  // Code for Arduino MEGA    
  DDRH = B11011111; //Set Digital PIN 8 as INPUT - All the others of port H as outputs
  PORTH = B00100000;
  DDRL = B00000100;//Set Digital PIN 47 as as OUTPUT - All the others of port L as inputs 
  PORTL = B00000000;
  noInterrupts();

}

void loop() {
  // put your main code here, to run repeatedly:
    if (PINH == B00000000){
    delayMicroseconds(33); //27us MODIFY THIS DELAY IF CHANGING THE SCAN PHASE
    PORTL = B00000100;
    delayMicroseconds(5);//5us
    PORTL = B00000000;
    delayMicroseconds(60);//60us
    PORTL = B00000100;
    delayMicroseconds(5);//5us
    PORTL = B00000000;
    
    }
}
