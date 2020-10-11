#include <fdeep/fdeep.hpp>

int main() {

	std::vector<float> inp;
	int cls_sz=6;
	float syn[]={1,7,2,2,1,1};

	inp.assign(syn, syn+cls_sz); 
	
    const auto model_0 = fdeep::load_model("DeepRL/fdeep_model.json");
        
    const fdeep::tensor inp_tnsr(fdeep::tensor_shape(static_cast<std::size_t>(cls_sz)),inp); //converting vector to tensor
    
    const auto result = model_0.predict({inp_tnsr});    
    std::cout << fdeep::show_tensors(result) << std::endl;
    
	const auto y=result[0];
    const std::vector<float> vec = y.to_vector();
    
    std::cout<<"vec[0], vec[1]: "<<vec[0]<<", "<<vec[1];
}
